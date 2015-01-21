function G = wanderlust(data, Options)
%    G = wanderlust( data, Options) % TODO structurize the return variable
%
%    Options                    : structure for specifying algorithm options which can contain
%                                   zero or more of the following fields:
%
%      [k]                      : truncate the neighbors to random k out of l 
%      [l]                      : size of neighborhood l>k closest points
%      [s]                      : index to the starting point in the data
%      [num_graphs]             : number of repeats - random selections of
%                                   k out of l.
%      [num_landmarks]          : number of waypoints\landmarks <- TODO
%                                 ensure logic if someone specifies predefined landmarks
%      [verbose]                : messages are printed to Matlab's sto
%      [metric]                 : string - distance metric for constructing 
%                                   the nearest neighbor graphs
%      [voting_scheme]          : How to weigh each point's contribution
%           'uniform'           -  
%           'exponential'       - 
%           'linear'            - 
%           'inverse-linear'    - 
%           TOOD return k_power -
%      [branch]                 : true\false to seek a branch
%      [band_sample]            : true\false landmarks are subsamples at equidistanced
%                                   bands away from the start point. using 
%                                   the points shortest path
%                                   distance over the graph.
%      [partial_order]          : an array of indices in data. the indices 
%                                point to landmarks and their order is forced 
%                                in the final output. 
%      [flock_landmarks]        : number of how many times to use a median
%      [plot_data]              : 2D data matrix to plot results to.
%                                  filter on landmarks
%      [snn]                    : shared nearest neighbor
%      [ann]                    : TODO - finalize - using Approx NN
%      [search_connected_components] :TODO add the option to cancel
%      [lnn] :DEBUG precomputed lnn
%      [landmarks] :DEBUG pre-chosen landmarks
%      
%
% run the wonderlust algorithm on data. generate num_graphs klNN graphs; starting point is s, choose num_landmarks
% random landmarks. weigh landmark-node pairs by distance to the power of power_k.
%
% (use distance as distance metric)
%
% (alternatively data can be lnn graph)
% 

% set up return structure
G.landmarks = [];
G.T = []; % traj
G.B = []; % branch

% Algorithm defaults
G.Opts.metric = 'cosine';
G.Opts.k = 8;
G.Opts.l = 30;
G.Opts.num_graphs = 25;
G.Opts.s = randsample(1:size(data, 1), 1);
G.Opts.num_landmarks = min(size(data, 1), 100);
G.Opts.verbose = true;
G.Opts.branch = false;
G.Opts.partial_order = [];
G.Opts.deblur = false;
G.Opts.snn = 0;
G.Opts.ann = false;
G.Opts.voting_scheme = 'linear';
G.Opts.band_sample = true;
G.Opts.flock_landmarks = 2;
G.Opts.search_connected_components = true;
G.Ops.plot_landmark_paths = false;
G.Opts.plot_data = [data(:,1) data(:,2)];
G.Opts.lnn = [];
G.Opts.landmarks = [];
G.Opts.disallow = [];


fn = fieldnames(Options);
for j=1:length(fn)
    name = fn{j};
    value = getfield(Options, name);
    
    if     strcmpi(name,'metric')           G.Opts.metric = value;
    elseif strcmpi(name,'k')                G.Opts.k = value;
    elseif strcmpi(name,'l')                G.Opts.l = value;
    elseif strcmpi(name,'num_graphs')       G.Opts.num_graphs = value;
    elseif strcmpi(name,'s')                G.Opts.s = value;
    elseif strcmpi(name,'num_landmarks')    G.Opts.num_landmarks = value;
    elseif strcmpi(name,'verbose')          G.Opts.verbose = value;
    elseif strcmpi(name,'branch')           G.Opts.branch = value;
    elseif strcmpi(name,'partial_order')   	G.Opts.partial_order = value;
    elseif strcmpi(name,'deblur')           G.Opts.deblur = value;
    elseif strcmpi(name,'snn')              G.Opts.snn = value;
    elseif strcmpi(name,'ann')              G.Opts.ann = value; %not supported yet
    elseif strcmpi(name,'voting_scheme')    G.Opts.voting_scheme = value; 
    elseif strcmpi(name,'band_sample')      G.Opts.band_sample = value; 
    elseif strcmpi(name,'flock_landmarks')  G.Opts.flock_landmarks = value; 
    elseif strcmpi(name,'plot_landmark_paths') G.Opts.plot_landmark_paths = value; 
    elseif strcmpi(name,'plot_data')        G.Opts.plot_data = value; 
    elseif strcmpi(name,'lnn')              G.Opts.lnn = value; 
    elseif strcmpi(name,'landmarks')        G.Opts.landmarks = value; 
    elseif strcmpi(name,'disallow')        G.Opts.disallow = value; 
    else   fprintf('Wanderlust.m: invalid option "%s" ignored.\n', name);
    end
end

% Build lNN graph
if issparse( data ) 
    if G.Opts.verbose 
        disp 'using prebuilt lNN graph';
    end
    lnn = data;
else
    if G.Opts.verbose 
        disp 'building lNN graph';
    end

    tic
    if false % tmp todo compute using the diffusion map code
        GraphDiffOpts = struct( ...
            'Normalization','smarkov', ...
            'Distance','euclidean', ...
            'Epsilon',1, ...
            'kNN', G.Opts.l, ...
            'kNNAutotune', 10, ... %'kEigenVecs', 6, ...
            'Symmetrization', 'W+Wt', ...
            'DontReturnDistInfo', 1 );

        GD = GraphDiffusion(data', 0, GraphDiffOpts);  
        lnn= GD.T;
    else
    if (isempty(G.Opts.lnn))
        lnn = parfor_spdists_knngraph( data, G.Opts.l,...
            'distance', G.Opts.metric,...
            'chunk_size', 500,... % TODO: parameterize and add opt for ppl without PC toolbox
            'verbose', G.Opts.verbose );
    else
        lnn = G.Opts.lnn;
    end
    end
    sprintf('lnn computed: %gs', toc);
    
    if (G.Opts.deblur)
        [i, j, s] = find(lnn);
        % flock each data point to its knn median
        for ith=1:numel(i)
            data(ith, :) = median(data(j(i==ith), :)); 
        end
        tic;
    	lnn = parfor_spdists_knngraph( data, l, 'distance', metric, 'chunk_size', 1000, 'SNN', true, 'verbose', true );
        sprintf('lnn re-computed after deblur: %gs', toc);
    end

    if ~isempty(G.Opts.disallow)

        [j, i, s] = find(lnn);
      
        % for each point 1..n
        nData = size(data,1);
        rem = cell(0, nData);
        parfor ci=1:nData
            
            % grab neighbors
            from = (ci-1)*G.Opts.l+1;
            to = ci*G.Opts.l;
            i_inds = from:to;
            i_group = G.Opts.disallow(ci);
            
            if i_group ~=0
            % for each neighbor
            for i_ind=i_inds
                i_neigh=j(i_ind);
                
                if G.Opts.disallow(i_neigh) ~=i_group || G.Opts.disallow(i_neigh)~=0 
                    
                    % add them to remove list
                    rem{ci} = [rem{ci} i_ind];
                end
            end
            end
        end

        rem = cell2mat(rem);
        
        % remove relevant indices
        i(rem) = [];
        j(rem) = [];
        s(rem) = [];
        lnn = sparse(j, i, s);

    end
    
    if (G.Opts.snn~=0)
        [j, i, s] = find(lnn);
        % observational note: i is sorted with l-1 apearences each index
        % use this irreliable observation to make sn faster
        
        nData = size(data,1);
        rem = cell(0, nData);

        % for each point 1..n
        parfor ci=1:size(data,1)
            
            % grab neighbors
            
            from = (ci-1)*G.Opts.l+1;
            to = ci*G.Opts.l;
            i_inds = from:to;
            i_neighs = j(i_inds);
            
            % for each neighbor
            for i_ind=i_inds
                i_neigh=j(i_ind);
                
                % grab the neighbor's neighbors
                from = (i_neigh-1)*G.Opts.l+1;
                to = i_neigh*G.Opts.l;
                j_neighs = j(from:to);
%                 j_neighs = j(i==i_neigh);
                
                % if the shared neighbors are not many enough
                if sum(ismember(i_neighs, j_neighs)) < G.Opts.snn
                    
                    % add them to remove list
                    rem{ci} = [rem{ci} i_ind];
                end
            end
        end
        
        rem = cell2mat(rem);

        % remove relevant indices
        i(rem) = [];
        j(rem) = [];
        s(rem) = [];
        lnn = sparse(j, i, s);
    end
end

% generate klNN graphs and iteratively refine a trajectory in each
G.klnn = {}
for graph_iter = 1:G.Opts.num_graphs
	if( G.Opts.verbose )
		fprintf( 1, 'iter #%d:\n', graph_iter );
	end

	iter_t = tic;

	% randomly generate a klNN graph
	if G.Opts.verbose
		fprintf( 1, 'entering knn graph: ' );
	end
   
    if (G.Opts.k~=G.Opts.l)
    	klnn = spdists_klnn( lnn, G.Opts.k, G.Opts.verbose );
        klnn = spdists_undirected( klnn ); % TODO consider removing - outliers?
    else
        klnn = lnn;
    end
 
	if( G.Opts.verbose )
		fprintf( 1, ' done (%3.2fs)\n', toc( iter_t ) );
		fprintf( 1, 'entering trajectory landmarks: ' );
	end

	% run traj. landmarks
	[ traj, dist, iter_l, RNK,paths_l2l ] = trajectory_landmarks( klnn,data, G);
    G.landmarks(graph_iter, :) = iter_l;
    G.traj(graph_iter) = {traj};
    G.dist(graph_iter) = {dist};
    G.klnn(graph_iter) = {klnn}; 
    
    if G.Opts.verbose
        fprintf( 1, ' done (%3.2fs)...\n', toc( iter_t ) );
    end

	% calculate weighed trajectory
	W = dist;
    if strcmpi(G.Opts.voting_scheme, 'uniform')
        W(:, :) = 1;
    elseif strcmpi(G.Opts.voting_scheme, 'exponential')
        sdv = mean ( std ( dist) )*3;
        W = W ./ repmat( sum( W ), size( W, 1 ), 1 );
        W = exp( -.5 * (W / sdv).^2);
    elseif strcmpi(G.Opts.voting_scheme, 'linear')
        W = W ./ repmat( sum( W ), size( W, 1 ), 1 );
        W = 1-W;
    end
    
    % The weghing matrix must be a stochastoc operator
    W_orig = W ./ repmat( sum( W ), size( W, 1 ), 1 );
    
    if (G.Opts.branch)
        cl = setdiff(unique(RNK)', RNK(G.Opts.s));
        W = W_orig;
        W(ismember(iter_l,find(RNK==cl(1))), RNK==cl(2)) = 0;
        W(ismember(iter_l,find(RNK==cl(2))), RNK==cl(1)) = 0;
        W = W ./ repmat( sum( W ), size( W, 1 ), 1 );
    end
    
    t( 1,:)  = traj(1,:);

	t( 2, : ) = sum( traj .* W );
	a = sum(W, 1);
    s = sum(a.*t(end, :));
    fprintf('s: %g\n', s);
    fprintf( 1, 'corr=');
    if any(isnan(t(2,:)))
            strPrompt = sprintf('there are nans in t, should I stop?');
            str_input = input(strPrompt,'s');
            if strcmpi(str_input, 'y') || strcmpi(str_input, 'yes') || strcmpi(str_input, 'yea')
                data(find(traj(2, :) == inf), :)
                return;
            end 
    end
    
	% iteratively realign trajectory (because landmarks moved)
	converged = 0; user_break = 0; realign_iter = 2;
	while( ~converged && ~user_break)
		realign_iter = realign_iter + 1;

		traj = dist;
		for idx = 1:size( dist, 1 )
			% find position of landmark in previous iteration
			idx_val = t( realign_iter - 1, iter_l( idx ) );
			% convert all cells before starting point to the negative
			before_indices = find( t( realign_iter - 1, : ) < idx_val );
			traj( idx, before_indices ) = -dist( idx, before_indices );
			% set zero to position of starting point
			traj( idx, : ) = traj( idx, : ) + idx_val;
        end

        if (G.Opts.branch)
            RNK = splittobranches(traj, t( end, : ),data, iter_l, dist,paths_l2l, G.Opts);
            cl = setdiff(unique(RNK)', RNK(G.Opts.s));
            if numel(cl) == 2
                W = W_orig;
                W(ismember(iter_l,find(RNK==cl(1))), RNK==cl(2)) = 0;
                W(ismember(iter_l,find(RNK==cl(2))), RNK==cl(1)) = 0;
                W = W ./ repmat( sum( W ), size( W, 1 ), 1 );
            else 
                fprintf( 'warning: no branch found');
            end
        end
        
		% calculate weighed trajectory
		t( realign_iter, : ) = sum( traj .* W );

        a = sum(W, 1);
%         s = sum(a.*t(end, :));
%         fprintf('s: %g\n', s);

		% check for convergence
        fprintf( 1, '%2.5f...', corr( t( realign_iter, : )', t( realign_iter - 1, : )' ));
		converged = corr( t( realign_iter, : )', t( realign_iter - 1, : )' ) > 0.9999;
        
        if (mod(realign_iter,40)==0)
            % !!!!temp!!!! break after too many alignments!
            user_break = true;
            
%             last_corr = corr( t( realign_iter, : )', t( realign_iter - 1, : )' );
%             strPrompt = sprintf('iter=%g, corr=%g. should I stop (y\n)?', realign_iter, last_corr);
%             str_input = input(strPrompt,'s');
%             if strcmpi(str_input, 'y') || strcmpi(str_input, 'yes') || strcmpi(str_input, 'yea')
%                 user_break = true;
%             end
        end
	end
	fprintf( 1, '\n%d realignment iterations, ', realign_iter );

	% save final trajectory for this graph    
    G.T(graph_iter, :) = t(realign_iter, :);

    if (G.Opts.branch)
        G.B(graph_iter, :) = RNK;
    end
    
    if( G.Opts.verbose )
		toc( iter_t );

		fprintf( 1, '\n' );
	end
end
end


function spdists = spdists_klnn( spdists, k, verbose )
	% spdists = spdists_klnn( spdists, k, verbose )
	%
	% given a lNN graph spdists, choose k neighbors randomly out of l for each node

	remove_edges = [];

	for idx = 1:length( spdists )
		% remove l-k neighbors at random
		neighs = find( spdists( :, idx ) );
		l = length( neighs ); % count number of neighbors
		remove_indices = neighs( randsample( length( neighs ), l - k ) );
		idx_remove_edges = sub2ind( size( spdists ), remove_indices, ones( l - k, 1 ) * idx );
		remove_edges = [ remove_edges; idx_remove_edges ];

		if( verbose )
			if( mod( idx, 50000 ) == 0 )
				fprintf( 1, '%3.2f%%', idx / length( spdists ) * 100 );
			elseif( mod( idx, 10000 ) == 0 )
				fprintf( 1, '.' );
			end
		end
	end

	spdists( remove_edges ) = 0;
    end

	function [ traj, dist, l, RNK,paths_l2l ] = trajectory_landmarks( spdists,data, G)
	% [ traj, dist, l ] = trajectory_landmarks( spdists, s, n, verbose )
	%
	% calculate the trajectory score of each point in spdists.
	%
	% s: list of indices of possible starting points. one of these points will be used to generate a reference
	% trajectory; the landmark shortest paths will be aligned to this reference.
	% n: list of landmark indices to use; or, alternatively, the number of landmarks to choose randomly from all
	% points.
	%
	% traj is a |n|x|spdists| matrix, where row i is the aligned shortest path from landmark i to each other point.
	% dist is a |n|x|spdists| matrix, where row i is the shortest path from landmark i to each other point. l is
	% the list of landmarks, l(1) is the starting point.

    RNK = zeros(size(data, 1), 1);
    n = G.Opts.num_landmarks;

    if( length( G.Opts.s ) > 1 )
		% if given a list of possible starting points, choose one. TODO move
		% to beginning of algorithm!!!!
		G.Opts.s = randsample( G.Opts.s, 1 );
	end

	if( length( n ) == 1 )
        [dists, paths, ~] = graphshortestpath( spdists, G.Opts.s,'METHOD','Dijkstra', 'directed', true);
        
        % if not given landmarks list, decide on random landmarks        
        if (G.Opts.band_sample)
            n_opts = [];
            num_jumps_arr = cellfun(@(x)numel(x), paths);
            max_jumps = max(num_jumps_arr);
            max_dist = max(dists);
            for prc = .998:-.10:.198
                band = find(dists>=(prc-.1)*max_dist & dists <=prc*max_dist); % & num_jumps_arr >= floor((prc-.05)*max_jumps) & num_jumps_arr <= ceil(prc*max_jumps));
                if length(band)> (n- 1 - length(G.Opts.partial_order))
                    n_opts = [n_opts randsample( band, n - 1 - length(G.Opts.partial_order), true )];
                end
            end
            n = randsample( n_opts, n - 1 - length(G.Opts.partial_order) );
        else
            n = randsample( 1:size(data,1), n - 1 - length(G.Opts.partial_order) );
        end
%         nw(end+1) = knnsearch(data, [-.2 -.2]);
%         nw(end+1) = knnsearch(data, [-.18 -.16]);
%         nw(end+1) = knnsearch(data, [.2 .5]);
%         nw(end+1) = knnsearch(data, [.1 .15]);
%         nw(end+1) = knnsearch(data, [.17 0]);
%         nw(end+1) = knnsearch(data, [-.1 .5]);
%         h=figure('Color',[1 1 1]);
%         [~, density, d1, d2] = kde2d(data, 1024);
%         contour(d1, d2, density, 128);
%         colormap(jet);
%         hold on;
%         plot(data(nw, 1),...
%             data(nw, 2),...
%             'Xm', 'markersize', 12,...
%             'linewidth',8);
%         plot(data(nw(end), 1),...
%             data(nw(end), 2),...
%             'Xm', 'markersize', 12,...
%             'linewidth',8);
        
        % flock landmarks 
        if (G.Opts.flock_landmarks > 0)
        for k=1:G.Opts.flock_landmarks
            [IDX, ~] = knnsearch(data, data(n, :), 'distance', G.Opts.metric, 'K', 20);     
            for i=1:numel(n)
                n(i) = knnsearch(data, median(data(IDX(i, :), :)), 'distance', G.Opts.metric); 
            end
        end
        end
    end

    partial_order = [G.Opts.s;G.Opts.partial_order(:)]; % partial_order includes start point
	l = [ partial_order; n(:) ]; % add extra landmarks if user specified
    
    % prune out weak edges
    prune = false;
    if (prune)
        res = 256;
        [band, density, x, y] = kde2d(data, res);

            % maybe remove edges jumping over low density regions?
            % BW1 = edge(density,'canny', .001);
            % imshow(BW1)
            % set(gca,'YDir','normal')

            %remove bad edges if the path to a landmark has a suspicious jump
        mapX = round( interp1(diag(x), 1:res, data(:, 1), 'linear', 'extrap') );
        mapY = round( interp1(diag(y), 1:res, data(:, 2), 'linear', 'extrap') );
        recalc = true;
        cou = 0;
        while recalc && cou < 30
            cou = cou+1;
            [dist( 1, : ), paths, ~] = graphshortestpath( spdists, s);%, 'directed', false );
            recalc = false;
            for pathidx=2:numel(l) %iterate over path to each landmark
                l_path = paths{l(pathidx)};
                path_jumpsX = abs(mapX(l_path(2:end))-mapX(l_path(1:end-1)));
                path_jumpsY = abs(mapY(l_path(2:end))-mapY(l_path(1:end-1)));
                path_jumps = path_jumpsX + path_jumpsY;
                bad_nodes = find(path_jumps > (mean(path_jumps) + 4*std(path_jumps)));
                if any(bad_nodes)
                    disp(sprintf('removing %g bad connections\n', numel(bad_nodes)));
                    spdists(sub2ind(size(spdists), l_path(bad_nodes), l_path(bad_nodes+1))) = 0;
                    spdists(sub2ind(size(spdists), l_path(bad_nodes+1), l_path(bad_nodes))) = 0;
                    recalc = true;
                end
            end
        end
    end    

	% calculate all shortest paths
    paths_l2l = cell(length(l));
    for li = 1:length( l )
        [dist( li, : ), paths, ~] = graphshortestpath( spdists, l( li ),'METHOD','Dijkstra');%, 'directed', false );
        paths_l2l(li) = {paths(l)};
        unreachable = (dist(li,:)==inf);
        while (any(unreachable) && G.Opts.search_connected_components)
            fprintf(['\n Warning: %g were unreachable. try increasing l'...
                'or k.Your data is possibly non continous, ie '...
                'has a completely separate cluster of points.'...
                'Wanderlust will roughly estimate their distance for now \n'],...
                sum(unreachable));
            % find closest unreachable point to reachable points.
            % connect it on the spdists. continue iteratively.
            unreachablei = find(unreachable);
            reachablei = find(~unreachable);
            while ~isempty(unreachablei)
                
                [idx, d] = knnsearch(data(unreachablei, :), data(reachablei, :));
                closest_reachable = d==min(d);
                
                %add connection to spdists
                spdists(reachablei(closest_reachable),...
                    unreachablei(idx(closest_reachable))) = min(d);
                spdists(unreachablei(idx(closest_reachable)),...
                    reachablei(closest_reachable)) = min(d);
                % move points from unreachable list to reachable
                reachablei(end+1:end+length(find(closest_reachable))) = ...
                    unreachablei(idx(closest_reachable));
                unreachablei(idx(closest_reachable)) = [];
                
                if sum(unreachable) > 100
                    break;
                end
            end
            [dist( li, : ), paths, ~] = graphshortestpath( spdists, l( li ),'METHOD','Dijkstra');%, 'directed', false );
            paths_l2l(li) = {paths(l)};
            unreachable = (dist(li,:)==inf);
        end
        
        if( G.Opts.verbose )
            fprintf( 1, '.' );
        end
    end
    
    
    % adjust paths according to partial order by redirecting
    nPartialOrder = length(partial_order);
    for radius = 1:nPartialOrder 
        for landmark_row = 1:nPartialOrder
            if (landmark_row + radius <= nPartialOrder)
                a = landmark_row;
                b = landmark_row + (radius-1);
                c = landmark_row + radius;
                dist(a, partial_order(c)) = dist(a, partial_order(b)) + dist(b, partial_order(c));
            end
            if (landmark_row - radius >= 1)
                a = landmark_row;
                b = landmark_row - (radius-1);
                c = landmark_row - radius;
                dist(a, partial_order(c)) = dist(a, partial_order(b)) + dist(b, partial_order(c));
            end
        end
    end

	% align to dist_1 - this for loop refers to partial order stuff
	traj = dist;
    for idx = 2:length(partial_order)
        [~, closest_landmark_row] = min(dist); %closest landmark will determine directionality
        traj(idx, closest_landmark_row < idx) = -dist(idx, closest_landmark_row < idx);
        traj( idx, : ) = traj( idx, : ) + dist( 1, l( idx ) );
    end
    
    % This is the actual align for regular wanderlust
    if length( l ) > length(partial_order)
        for idx = length(partial_order)+1:length( l )
            % find position of landmark in dist_1
            idx_val = dist( 1, l( idx ) );
            % convert all cells before starting point to the negative
            before_indices = find( dist( 1, : ) < idx_val );
            traj( idx, before_indices ) = -dist( idx, before_indices );
            % set zero to position of starting point
            traj( idx, : ) = traj( idx, : ) + idx_val;
        end
    end
    if (G.Opts.plot_landmark_paths)
        plot_landmark_paths(G.Opts.plot_data, paths_l2l, l);
    end
    if (G.Opts.branch)
        [RNK, bp] = splittobranches(traj, traj(1, :), data, l, dist, paths_l2l, G.Opts);
    end
end

function [RNK, pb] = splittobranches(trajs, t, data, landmarks, dist, paths_l2l, Opts)
    
    proposed = repmat(t(landmarks), size(trajs, 1), 1);
    reported = trajs(1:length(landmarks), landmarks);

    % square matrix of the difference of perspectives landmark to landmark
    diffdists = abs(reported - proposed);
    diffdists = .5*(diffdists'+diffdists);
    c = segmentlikemichealjordanwould(diffdists);
    
    % show wine glass with clusterization
    [evec2, ~] = eig(diffdists);
    evec2 = evec2(:, 2);
    [~, idx] = sort(evec2);
%     figure('Color',[1 1 1]);
%     scatter(evec2(idx),t(landmarks(idx)), ones(size(evec2))*50, c(idx));
%     figure('Color',[1 1 1]);
%     subplot(1,2,1);
%     imagesc(diffdists(idx,idx));
%     subplot(1,2,2);
%     [~, idx_time] = sort(t(landmarks));
%     imagesc(diffdists(idx_time,idx_time));
    
    % to pinpoint branch - look into the min (traj) of the paths from one cluster to the
    % other
    c_branch = setdiff(unique(c)', c(1)); % the branch indices
    paths_branch_a = paths_l2l(c==c_branch(1));
    paths_branch_b = paths_l2l(c==c_branch(2));
    fork_p = [];
    for i=1:numel(paths_branch_a)
        paths_branch_a_to_b = paths_branch_a{i}(c==c_branch(2));
        for j=1:numel(paths_branch_a_to_b)
            fork_p(end+1) = min(t(paths_branch_a_to_b{j}));
        end
    end
    
    for i=1:numel(paths_branch_b)
        paths_branch_b_to_a = paths_branch_b{i}(c==c_branch(1));
        for j=1:numel(paths_branch_b_to_a)
            fork_p(end+1) = min(t(paths_branch_b_to_a{j}));
        end
    end
    
    % reassign to clusters based on branch point
    pb = prctile(fork_p, 10);
    c_new = c;
    c_new(t(landmarks)' <= pb) = c(1);
    c_new(evec2<0 & t(landmarks)' >= pb) = c_branch(1);
    c_new(evec2>0 & t(landmarks)' >= pb) = c_branch(2);
    
    if (numel(unique(c_new))<3)
        % we have an empty branch
        % probably some outlier population is taking over the v2 signal
        fprintf('problems:');
        figure('Color',[1 1 1]);
        subplot(1,2,1);
        
        scatter(evec2(idx),t(landmarks(idx)), ones(size(evec2))*50, c(idx));
        title(sprintf('Problems: BP=%g', pb));
        
        subplot(1,2,2);
        scatter(Opts.plot_data(:,1),Opts.plot_data(:,2),...
            ones(size(data,1),1)*10, ones(size(data,1),1)); 
        hold on;
        scatter(Opts.plot_data(landmarks,1),Opts.plot_data(landmarks,2),...
            ones(numel(landmarks),1)*50, (c+1)*5);
        
    end
        
%     figure('Color',[1 1 1]);
%     scatter(evec2,t(landmarks), ones(size(evec2))*50, c_new);
    
    % what do we get if we michal jordan the wine glass instead?
%     figure('Color',[1 1 1]);
%     c_glass = segmentlikemichealjordanwould([(evec2-mean(evec2))/std(evec2) (trajs(1,landmarks)-mean(trajs(1,landmarks)))'/(std(trajs(1,landmarks))^2)]);
%     scatter(evec2(idx),trajs(1,landmarks(idx)), ones(size(evec2))*50, c_glass(idx));
    
    % what about denoising diffdists
%     [L,S] = ADMMRPCA(diffdists(idx,idx)./max(max(diffdists)), 0.1,1); 
    
    % for each datapoint find closest landmark
    [~,I] = min(abs(dist(1:numel(landmarks), :)));
    RNK = c_new(I);
end

function c=segmentlikemichealjordanwould(data)
sigma = 10000;
sigma = prctile(data(:), 97);

% form affinity matrix with gaussian eucliean distance
A = exp((-.5*(1/sigma^2)).*pdist2(data, data).^2);

% compute the laplacian
D = diag(sum(A,2).^(-.5));
L = D*A*D;

% find x1,x2,x3 the 3 largest eigenvector of L (chosen to be orthogonal in
% case of repeated?)
[evec, eval] = eig(L);
X = evec(:, [1 2 3]);
Y = X./repmat(sqrt(sum(X.^2,2)), 1, size(X, 2));
c = kmeans(Y, 3);

end

function plot_landmark_paths(data, paths, l)
    figure('Color',[1 1 1]);
    scatter (data(:, 1), data(:, 2));
    hold on;
    plot(data(l(1), 1), data(l(1), 2), 'Xr');
    scatter(data(l, 1), data(l, 2), 20*ones(numel(l), 1), 'or');
    for p=1:numel(paths)
        pathsp = paths{p};
        for q=1:numel(pathsp)
            pathq=pathsp{q};
            plot(data(pathq, 1), data(pathq, 2), 'k-');
        end
    end
end