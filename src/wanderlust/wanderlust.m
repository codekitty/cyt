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
%      [search_connected_components] : search for connected components if
%                                   graph is discontious
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
G.Opts.plot_landmark_paths = false;
G.Opts.plot_data = [data(:,1) data(:,2)];
G.Opts.lnn = [];
G.Opts.landmarks = [];
G.Opts.disallow = [];
G.Opts.cell_clusters = [];
G.Opts.end_clusters = [];
G.Opts.plot_debug_branch = true;

rng('shuffle');

fn = fieldnames(Options);
for j=1:length(fn)
    name = fn{j};
    
    G.Opts.(deblank(name)) = Options.(deblank(name));
    
%     value = getfield(Options, name);
%     if     strcmpi(name,'metric'),          G.Opts.metric = value;
%     elseif strcmpi(name,'k'),               G.Opts.k = value;
%     elseif strcmpi(name,'l'),               G.Opts.l = value;
%     elseif strcmpi(name,'num_graphs'),      G.Opts.num_graphs = value;
%     elseif strcmpi(name,'s'),               G.Opts.s = value;
%     elseif strcmpi(name,'num_landmarks'),   G.Opts.num_landmarks = value;
%     elseif strcmpi(name,'verbose'),         G.Opts.verbose = value;
%     elseif strcmpi(name,'branch'),          G.Opts.branch = value;
%     elseif strcmpi(name,'partial_order'),  	G.Opts.partial_order = value;
%     elseif strcmpi(name,'deblur'),          G.Opts.deblur = value;
%     elseif strcmpi(name,'snn'),             G.Opts.snn = value;
%     elseif strcmpi(name,'ann'),             G.Opts.ann = value; %not supported yet
%     elseif strcmpi(name,'voting_scheme'),   G.Opts.voting_scheme = value; 
%     elseif strcmpi(name,'band_sample'),     G.Opts.band_sample = value; 
%     elseif strcmpi(name,'flock_landmarks'), G.Opts.flock_landmarks = value; 
%     elseif strcmpi(name,'plot_landmark_paths'), G.Opts.plot_landmark_paths = value; 
%     elseif strcmpi(name,'plot_data'),       G.Opts.plot_data = value; 
%     elseif strcmpi(name,'lnn'),             G.Opts.lnn = value; 
%     elseif strcmpi(name,'landmarks'),       G.Opts.landmarks = value; 
%     elseif strcmpi(name,'disallow'),        G.Opts.disallow = value; 
%     elseif strcmpi(name,'search_connected_components'), G.Opts.search_connected_components = value;
%     elseif strcmpi(name,'cell_clusters'),   G.Opts.cell_clusters = value;
%     elseif strcmpi(name,'end_clusters'),    G.Opts.end_clusters = value;
%     else   fprintf('Wanderlust.m: invalid option "%s" ignored.\n', name);
%     end
end

G.Opts.plot_debug_branch = false;

G.Opts

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

    tic;
    if false % tmp todo compute using the diffusion map code
        GraphDiffOpts = struct( ...
            'Normalization','smarkov', ...
            'Distance',G.Opts.metric, ...
            'Epsilon',1, ...
            'kNN', G.Opts.l, ...
            'kNNAutotune', 10, ... %'kEigenVecs', 6, ...
            'Symmetrization', 'W+Wt', ...
            'DontReturnDistInfo', 1 );

        GD = GraphDiffusion(data', 0, GraphDiffOpts);  
        lnn= GD.T;
    else
        
        if (isempty(G.Opts.lnn))

            % hash map for data+metric+l -> lnn
            hashMat = DataHash(data); 
            hashVal = DataHash(string2hash([G.Opts.metric hashMat]) + G.Opts.l);

            % Check for cached
            [curr_path, ~, ~] = fileparts(mfilename('fullpath'));
            curr_path = [curr_path filesep];
            cachefilename = [curr_path 'cacheknnResults.mat'];

            %checking for old lnn results for the same data
            fileCheck = exist (cachefilename, 'file');

            % if no file
            if (fileCheck==0) 
                % create new hash map
                mapMat = containers.Map(); 
            else
                % loading the old lnn results
                file= load(cachefilename); 
                mapMat=file.mapMat;
            end

            % check for key in the hash map
            check = isKey(mapMat,hashVal); 

            % if lnn found in cache
            if (check==1) 
                % no need to run lnn again->returning old result
                value=values(mapMat,{hashVal});
                lnn=value{1};
                fprintf('lnn loaded from cache: %gs\n', toc);
            else
                lnn = parfor_spdists_knngraph( data, G.Opts.l,...
                    'distance', G.Opts.metric,...
                    'chunk_size', 1000,... % TODO: parameterize and add opt for ppl without PC toolbox
                    'verbose', G.Opts.verbose );
                fprintf('lnn computed: %gs\n', toc);
                
                % while the hash map is too big removing the first element
                while length(mapMat)>5 
                    lstKeys= keys(mapMat); 
                    remove(mapMat,lstKeys(1));
                end

                % adding the name and lnn result to hashmap
                mapMat(hashVal)=lnn; 

                % saving into file
                save(cachefilename,'mapMat');
            end
        else
            lnn = G.Opts.lnn;
        end
    end
    
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
        rem = cell(1, nData);

        tic;
        % for each point 1..n
        parfor ci=1:nData
            
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
        sprintf('shared nearest neighbors computed: %gs', toc);

        rem = cell2mat(rem);

        % remove relevant indices
        i(rem) = [];
        j(rem) = [];
        s(rem) = [];
        lnn = sparse(j, i, s);
    end
end

% generate klNN graphs and iteratively refine a trajectory in each
G.klnn = {};
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
    else
        klnn = lnn;
    end
    
    G.Opts.exclude_points = [];
    G.Opts.minimum_incoming = 0;
    if G.Opts.minimum_incoming > 0 
       [i,j, s] = find(klnn);
        tabs = tabulate(i);
        inds = tabs((tabs(:,2)<G.Opts.minimum_incoming), 1);
        i_new = i(~ismember(j, inds));
        j_new = j(~ismember(j, inds));
        s_new = s(~ismember(j, inds));
        klnn = sparse(i_new, j_new, s_new);
        
        if (length(G.Opts.num_landmarks) > 1) 
            G.Opts.num_landmarks = setdiff(G.Opts.num_landmarks, inds);
        end
        
        G.Opts.exclude_points = inds;
        % TODO if G.Opts.minimum_incoming>1 then we need to remove also the
        % outgoing edges
    end
    
    klnn = spdists_undirected( klnn ); % TODO consider removing - outliers?
    G.lnn = klnn;

 
	if( G.Opts.verbose )
		fprintf( 1, ' done (%3.2fs)\n', toc( iter_t ) );
		fprintf( 1, 'entering trajectory landmarks: ' );
	end

	% run traj. landmarks
	[ traj, dist, iter_l, RNK,paths_l2l, diffdists,Y] = trajectory_landmarks( klnn,data, G);
    
    % save output variables
    G.landmarks(graph_iter, :) = iter_l;
    G.traj(graph_iter) = {traj};
    G.dist(graph_iter) = {dist};
    G.klnn(graph_iter) = {klnn}; 
    
    if G.Opts.verbose
        fprintf( 1, ' done (%3.2fs)...\n', toc( iter_t ) );
    end

	% calculate weighed trajectory
    if strcmpi(G.Opts.voting_scheme, 'uniform')
        W_full(:, :) = 1;
    elseif strcmpi(G.Opts.voting_scheme, 'exponential')
        sdv = mean ( std ( dist) )*3;
        W_full = exp( -.5 * (dist / sdv).^2);
    elseif strcmpi(G.Opts.voting_scheme, 'linear')
        W_full = repmat(max(dist), size( dist, 1 ), 1) - dist;
        if ~isempty(G.Opts.exclude_points)
            W_full(:, G.Opts.exclude_points) = 1;
        end
    end
    
    % The weghing matrix must be a column stochastic operator
    W_full = W_full ./ repmat( sum( W_full ), size( W_full, 1 ), 1 );
        
    if (G.Opts.branch)
        W = muteCrossBranchVoting(W_full, RNK, RNK(G.Opts.s), iter_l, Y);
    else
        W = W_full;
    end
    
    % save initial solution - start point's shortest path distances
    t=[];
    t( 1,:)  = traj(1,:);
	t( end+1, : ) = sum( traj .* W );
    
	% iteratively realign trajectory (because landmarks moved)
	converged = 0; user_break = 0; realign_iter = 2;
    G.v = [];
	while  ~converged && ~user_break
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
%             [RNK, bp, diffdists, Y] = splittobranches(traj, traj(1, : ),data, iter_l, dist,paths_l2l, G.Opts);
            W = muteCrossBranchVoting(W_full, RNK, RNK(G.Opts.s), iter_l,Y);
        end
       
        G.v(end+1) = compute_energy(t(end, iter_l), dist(:, iter_l));
        
		% calculate weighed trajectory
		t( realign_iter, : ) = sum( traj .* W );

		% check for convergence
        fpoint_corr = corr( t( realign_iter, : )', t( realign_iter - 1, : )' );
        fprintf( 1, '%2.5f...', fpoint_corr);
		converged = fpoint_corr > 0.9999;
        
        if (mod(realign_iter,16)==0)
            % break after too many alignments - something is wrong
            user_break = true;
            fprintf('\nWarning: Force exit after %g iterations\n', realign_iter);
        end
	end
    
    G.v(end+1) = compute_energy(t(end, iter_l), dist(:, iter_l));
    figure('Color',[1 1 1]);
    plot(1:numel(G.v), G.v, '-b');
%     plot_iterations(G.Opts.plot_data, t);
    
	fprintf( 1, '\n%d realignment iterations, ', realign_iter-1 );

	% save final trajectory for this graph    
    G.T(graph_iter, :) = t(realign_iter, :);
    
    if ~isempty(G.Opts.exclude_points)
        nearest_landmarks = knnsearch(data(iter_l, :), data(G.Opts.exclude_points, :));        
        G.T(graph_iter, G.Opts.exclude_points) = G.T(graph_iter, iter_l(nearest_landmarks));
    end
    
    if (G.Opts.branch)
        % Recalculate branches post reassignments
%         [RNK, bp, diffdists, Y] = splittobranches(traj, traj(1, : ), data, iter_l, ...
%             dist,paths_l2l, G.Opts);
bp=0;
        G.B(graph_iter, :) = RNK;
        G.diffdists = diffdists;
        G.bp(graph_iter) = bp;
        G.Y(graph_iter, :) = Y;
    else
        G.B = G.T; % branch
        G.bp(graph_iter) = 0;
    end
    
	if( G.Opts.verbose )
		toc( iter_t );

		fprintf( 1, '\n' );
	end
end
end


% spdists = spdists_klnn( spdists, k, verbose )
%
% given a lNN graph spdists, choose k neighbors randomly out of l for each node
% consider re-writing this using find on fll spdists and sparse to recreate
function spdists = spdists_klnn( spdists, k, verbose )

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
function [ traj, dist, l, RNK,paths_l2l, diffdists,Y ] = trajectory_landmarks( spdists,data, G)

    Y = [];
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
        n_opts = 1:size(data,1);
        if (G.Opts.band_sample)
            n_opts = [];
            window_size = .1;
            max_dist = max(dists);
            for prc = .998:-window_size:.08
                band = find(dists>=(prc-window_size)*max_dist & dists <=prc*max_dist); % & num_jumps_arr >= floor((prc-.05)*max_jumps) & num_jumps_arr <= ceil(prc*max_jumps));
                n_opts = [n_opts randsample( band, min(length(band), n - 1 - length(G.Opts.partial_order)), false )];
            end
        end
        n = randsample( n_opts, n - 1 - length(G.Opts.partial_order) );
        
        % if branch - add 'tailk' landmarks from the tail (30 prc) of the data
        if G.Opts.branch
            tailk=20;
            [dists, ~, ~] = graphshortestpath( spdists, G.Opts.s,'METHOD','Dijkstra', 'directed', true);
            tailband = find( dists>=(prctile(dists, 90)) );
            tailk = min([length(tailband), tailk, floor(length( n )/2)]); 
            n(randsample( 2:length( n ), tailk)) = randsample( tailband, tailk);
        end
        
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

    diffdists = zeros(length(n), length(n));

    partial_order = [G.Opts.s;G.Opts.partial_order(:)]; % partial_order includes start point
	l = [ partial_order; n(:) ]; % add extra landmarks if user specified
    
	if true
    % calculate all shortest paths
    paths_l2l = cell(length(l),1);
    for li = 1:length( l )
        [dist( li, : ), paths, ~] = graphshortestpath( spdists, l( li ),'METHOD','Dijkstra', 'directed', false );
        if sum(cellfun(@(x)isempty(x), paths(l))) 
            fprintf('\nWarning: found empty path');
        end
        paths_l2l(li) = {paths(l)};
        unreachable = (dist(li,:)==inf);
        unreachable(G.Opts.exclude_points) = 0;

        while (any(unreachable) && G.Opts.search_connected_components)
            fprintf(['\n Warning: %g were unreachable. try increasing l'...
                'or k.Your data is possibly non continous, ie '...
                'has a completely separate cluster of points.'...
                'Wanderlust will roughly estimate their distance for now \n'],...
                sum(unreachable));
            if (G.Opts.plot_debug_branch)
                figure('Color',[1 1 1]);
                scatter(G.Opts.plot_data(:,1), G.Opts.plot_data(:,2), 150,'.b');
                hold on;
                scatter(G.Opts.plot_data(l(li),1), G.Opts.plot_data(l(li),2), 150,'.g');
                scatter(G.Opts.plot_data(unreachable,1), G.Opts.plot_data(unreachable,2), 150,'.r');
                title('Unreachable in red (from gre');
            end
            % find closest unreachable point to reachable points.
            % connect it on the spdists. continue iteratively.
            unreachablei = find(unreachable);
            reachablei = find(~unreachable);
            cou = 0;
            while ~isempty(unreachablei)
                cou = cou+1;
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
                
                if ~mod(cou, 10)
                    break;
                end
            end
            [dist( li, : ), paths, ~] = graphshortestpath( spdists, l( li ),'METHOD','Dijkstra', 'directed', false );
            paths_l2l(li) = {paths(l)};
            unreachable = (dist(li,:)==inf);
        end
        
        if( G.Opts.verbose )
            fprintf( 1, '.' );
        end
    end
    if ~isempty(G.Opts.exclude_points)
        dist(:, G.Opts.exclude_points) = mean(mean(dist~=inf));
    end
    else

    kev = 15;
    GraphDiffOpts = struct( 'Normalization','smarkov', ...
                            'Epsilon',1, ...
                            'kNN', 15, ...
                            'kEigenVecs', kev, ...
                            'Symmetrization', 'W+Wt'); ...

    GD = GraphDiffusion(data', 0, GraphDiffOpts);
    paths_l2l = cell(length( l ),1);
    for li = 1:length( l )
        dist(li, :) = 100*sqrt(sum(bsxfun(@minus, GD.EigenVecs(:, 2:kev), GD.EigenVecs(l(li), 2:kev)).^2, 2));
        paths={};

       [~, paths, ~] = graphshortestpath( spdists, l( li ),'METHOD','Dijkstra', 'directed', false );
       paths_l2l{li} = paths(l);
        if( G.Opts.verbose )
            fprintf( 1, '.' );
        end
    end
    end
    
    if any(any(dist==inf))
        dist(dist==inf) = max(max(dist~=inf));
        if (G.Opts.verbose)
            fprintf('\nwarning: some points remained unreachable (dist==inf)');
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
%         plot_landmark_paths(G.Opts.plot_data, paths_l2l, l);
    end
    if (G.Opts.branch)
        [RNK, bp, diffdists, Y] = splittobranches(traj, traj(1, :), data, l, dist, paths_l2l, G.Opts);
    end
end

function [RNK, pb, diffdists, Y] = splittobranches(trajs, t, data, landmarks, dist, paths_l2l, Opts)
    
    proposed = repmat(t(landmarks), size(trajs, 1), 1);
    reported = trajs(1:length(landmarks), landmarks);

    % Extract landmark clusters if specified
    if (length(Opts.cell_clusters) > 0)
        landmark_clusters = Opts.cell_clusters(landmarks);
    else
        landmark_clusters = [];
    end


    % square matrix of the difference of perspectives landmark to landmark
    diffdists = abs(reported - proposed);
    diffdists = (diffdists'+diffdists)/2;
    
%     %normalize by dist?
%     diffdists = diffdists.*repmat(dist(1,landmarks)', 1, numel(landmarks));

%     c = segmentlikemichealjordanwould(diffdists);
%     c = Opts.end_clusters(landmarks);
    
    % get second eigen vector of diifdists
    [EigenVecs, EigenVals] = GetEigs(diffdists, 2, [],struct('TakeDiagEigenVals',1));
    evec2 = EigenVecs(:, 2);
    [~, idx] = sort(evec2);
    
    % assign last positive 5 and last negative 5
    t_l = t(landmarks);
    b1_suspect = find(evec2<0); % suspects for branch 1
    b2_suspect = find(evec2>0); % suspects for branch 2
    
    [~, b1_sorted_inds] = sort(t_l(b1_suspect), 'descend');
    [~, b2_sorted_inds] = sort(t_l(b2_suspect), 'descend');
    
    c = ones(1, numel(landmarks));
    c(b1_suspect(b1_sorted_inds(1:min(5, numel(b1_sorted_inds))))) = 2;
    c(b2_suspect(b2_sorted_inds(1:min(5, numel(b2_sorted_inds))))) = 3;
    
	trunk = c(1);

    c_branch = setdiff(unique(c)', trunk); % the branch indices
    brancha = find(c==c_branch(1));
    branchb = find(c==c_branch(2));
    paths_branch_a = paths_l2l(brancha);
    paths_branch_b = paths_l2l(branchb);
    fork_p = [];
    for i=1:numel(paths_branch_a)
        paths_branch_a_to_b = paths_branch_a{i}(branchb);
        for j=1:numel(paths_branch_a_to_b)
            if isempty(paths_branch_a_to_b{j})
                fprintf('no path from l:%g to l:%g', brancha(i), branchb(j));
            else
                fork_p(end+1) = min(t(paths_branch_a_to_b{j}));
            end
        end
    end
    
    for i=1:numel(paths_branch_b)
        paths_branch_b_to_a = paths_branch_b{i}(brancha);
        for j=1:numel(paths_branch_b_to_a)
           if isempty(paths_branch_b_to_a{j})
                fprintf('no path from l:%g to l:%g', branchb(i), brancha(j));
           else
                fork_p(end+1) = min(t(paths_branch_b_to_a{j}));
           end
        end
    end
    
    % reassign to clusters based on branch point
    pb = prctile(fork_p, 50);
    c_new = c;
    [~,I] = min(abs(dist(1:numel(landmarks), :)));
    RNK = c_new(I);
    c_new(t(landmarks)' <= pb) = c(1);
    c_new(evec2<0 & t(landmarks)' >= pb) = c_branch(1);
    c_new(evec2>0 & t(landmarks)' >= pb) = c_branch(2);

    e2 = evec2;
    timev = t(landmarks);
    Y = e2(I);
    
    % Compute affinity matrix over landmark distances
    n = size(dist,2); % num points
    regular_points = 1:n;%setdiff(1:n, landmarks);
    sigma = .25*std(dist(:));
    Aff = exp((-.5*(1/sigma^2)).*dist(:, regular_points).^2);
    
    % make aff matrix a stochastic operator
    Stoch=columndiv(Aff, sum(Aff));
    Y = zeros(1,n);
    Y(landmarks) = e2;
    Y(regular_points)=(Stoch'*e2(:)).*(t(regular_points)'.^.7);
    
    if (Opts.plot_landmark_paths && (Opts.plot_debug_branch || numel(unique(c_new))<3))

        fh = figure('Color',[1 1 1]);
        
        subplot(2,2,1);       
        scatter(Y,t, 150, '.');
        title('Projection');

        subplot(2,2,2);       
        scatter(evec2(idx),t(landmarks(idx)), ones(size(evec2))*150, c_new(idx), '.');
        title(sprintf('BP ident @ tau=%g', pb));

        subplot(2,2,3);
        scatter(Opts.plot_data(:,1),Opts.plot_data(:,2),...
            ones(size(data,1),1)*30, '.b'); 
        hold on;
        
        scatter(Opts.plot_data(landmarks(c==1),1),Opts.plot_data(landmarks(c==1),2),...
            ones(numel(landmarks(c==1)),1)*50, 'ok');
        scatter(Opts.plot_data(landmarks(c==2),1),Opts.plot_data(landmarks(c==2),2),...
            ones(numel(landmarks(c==2)),1)*50, 'or');
        scatter(Opts.plot_data(landmarks(c==3),1),Opts.plot_data(landmarks(c==3),2),...
            ones(numel(landmarks(c==3)),1)*50, 'og');
        
        subplot(2,2,4);
        scatter(Opts.plot_data(:,1),Opts.plot_data(:,2),...
            ones(size(data,1),1)*30, '.b'); 
        hold on;

        scatter(Opts.plot_data(landmarks(c_new==1),1),Opts.plot_data(landmarks(c_new==1),2),...
            ones(numel(landmarks(c_new==1)),1)*50, 'ok');
        scatter(Opts.plot_data(landmarks(c_new==2),1),Opts.plot_data(landmarks(c_new==2),2),...
            ones(numel(landmarks(c_new==2)),1)*50, 'or');
        scatter(Opts.plot_data(landmarks(c_new==3),1),Opts.plot_data(landmarks(c_new==3),2),...
            ones(numel(landmarks(c_new==3)),1)*50, 'og');
       
        print(f, '-dpng', sprintf('%giteration.png', timestamp), '-r100');
        
        % show Q sorted by second eig vector, 
        figure('Color',[1 1 1]);
        
        subplot(1,2,1);
        imagesc(diffdists(idx,idx));
        set(gca,'xtick',[],'ytick',[]);
        drawnow;
       	ax = gca;
        ax.YTickLabel = cellfun(@num2str, num2cell(c(idx)), 'UniformOutput', false);
        title('sorted by second eigen vector');
        colormap jet

        % show Q sorted by tau, 
        subplot(1,2,2);
        [~, idx_time] = sort(t(landmarks));
        imagesc(diffdists(idx_time,idx_time));
        set(gca,'xtick',[],'ytick',[]);
        drawnow;
       	ax = gca;
        ax.YTickLabel = cellfun(@num2str, num2cell(c(idx_time)), 'UniformOutput', false);
        title('sorted by tau');
        colormap jet
        drawnow;
    end
        
%     figure('Color',[1 1 1]);
%     scatter(evec2,t(landmarks), ones(size(evec2))*50, c_new);
    
    % what do we get if we michael jordan the wine glass instead?
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
    sigma = prctile(data(:), 97);

    % form affinity matrix with gaussian eucliean distance
    A = exp((-.5*(1/sigma^2)).*pdist2(data, data).^2);

    % compute the laplacian
    D = diag(sum(A,2).^(-.5));
    L = D*A*D;

    % Kmeans on normalized eigen vectors
    [evec, ~] = GetEigs(L, 3, [],struct('TakeDiagEigenVals',1));

    X = evec(:, 1:3);
    Y = X./repmat(sqrt(sum(X.^2,2)), 1, size(X, 2));
	c = kmeans(Y, 3);
end

function plot_landmark_paths(data, paths, l)
    figure('Color',[1 1 1]);
    nData = size(data, 1);
    scatter (data(:, 1), data(:, 2), 2*ones(1, nData), '.b');
    hold on;
    plot(data(l(1), 1), data(l(1), 2),20, 'Xr');
    scatter(data(l, 1), data(l, 2), 20*ones(numel(l), 1), 'or');
    for p=1:numel(paths)
        pathsp = paths{p};
        for q=1:numel(pathsp)
            pathq=pathsp{q};
            plot(data(pathq, 1), data(pathq, 2), 'k-');
        end
    end
    drawnow;
end

function W=muteCrossBranchVoting(W, RNK, trunk_id, landmarks, Y)

    % range between -1 and 1
    Y_scale = Y-median(Y(landmarks));
    Y_scale(Y_scale<0) = Y_scale(Y_scale<0)./max(abs((Y_scale(Y_scale<0))));
    Y_scale(Y_scale>0) = Y_scale(Y_scale>0)./max(Y_scale(Y_scale>0));
    Y_pos = abs(Y_scale)';

    W_test= W;
    b=std(Y_scale);
    for i=1:numel(Y)
        crossb = sign(Y_scale(i))~=sign(Y_scale(landmarks));
        W_test(crossb, i) = W(crossb,i).*...
                                max(exp(-.5*(Y_pos(i).^2)./b),...
                                    exp(-.5*(Y_pos(landmarks(crossb)).^2)./b));
    end   
    W = columndiv(W_test, sum( W_test ));
    
%     i = 47649 
%     [W(:,i) W_test(:, i) Y_scale(landmarks)' Y_scale(i)*ones(length(landmarks),1)]

%     % grab branch cluster labels
%     branch_ids = setdiff(unique(RNK)', trunk_id);
% 
%     % if we have indeed 2 branches
%     if numel(branch_ids) == 2
% 
%         % mute voting weight between one branch to the other
%         W(ismember(landmarks,find(RNK==branch_ids(1))), RNK==branch_ids(2)) = 0;
%         W(ismember(landmarks,find(RNK==branch_ids(2))), RNK==branch_ids(1)) = 0;
% 
%         % make W column stochastic (weighted average per landmark)
%         W = W ./ repmat( sum( W ), size( W, 1 ), 1 );
% 
%     % otherwise, print warning
%     else 
%         fprintf( 'warning: no branch found');
%     end
end

function v=compute_energy(tau, D)
        n = numel(tau);
        sdv = mean ( std ( D) )*3;
        E = exp( -.5 * (D / sdv).^2);
        L = E-diag(sum(E, 1));
        
        % compute error
        v = 0;
        for i=1:n
            for j=1:n
                v = v - L(i,j)*((tau(i)-tau(j))^2);
            end
        end
        
        for i=1:n
            for j=1:n
                v = v - 2 * L(i,j)*(abs(tau(i)-tau(j)))*D(i,j);
%                 v = v + E(i,j)*(abs(tau(i)-tau(j)))*D(i,j);
            end
        end
        
        % assert sum is constant
        lambda = sum(E, 1);
        s = sum(lambda.*tau);        
        
        fprintf('\nV = %g\n S = %g\n',v, s);
end

function plot_iterations(data, t)
    for iter=1:size(t, 1)
        % new white background figure
        figure('Color',[1 1 1]);
        
        % scatter the data with solution on top
        scatter(data(:,1),data(:,2), 35, t(iter,:), '.'); 
        
        % make it look nice
        title(sprintf('iteration: %g', iter));
        colormap jet;
        box on;
        set(gca,'xtick',[],'ytick',[]);
    end
end