function plot_as_function(x, Y, varargin)
% plots Y (NxM) as a function of vector x (N length).  
% an average of points in Y is computed for 'num_locs' along vector x.
%
% 'labels': a string cell array to enter for the legend
%
% specify 'avg_type': default 'gaussian_var'
%
% 'sliding': a fixed averaged is computed within a sliding window. takes
% into account points within two windows to each side. so let's say 100
% locations would mean each location is averaging points from two windows to
% the right and to the left
%
% 'linear': same as sliding the each point's contibution dimishes linearly
% 
% 'squared': same as sliding the each point's contibution dimishes
% polynomialy
%
% 'gaussian': all points contibute to the weighted average but beyond 2
% standard deviations the contibution is insignificant.
%
% 'gaussian_var': same as gaussian but the standard deviation is variable
% as a function of the normalized density of the points along x. so in
% saprser locations on x, we use a larger standard deviation to let farther
% points contribute.
%
%
% 'show_error': true or false, shows error bars along the curve, default
% 'false'
% 
% Michelle Tadmor, Columbia University, 2013-2015

Markers={'--','-', ':'};

clear persistent_ksdensity;
legend_flag = false;
avg_type = 'gaussian_var';
num_locs = 100;
show_error = false;
normalize = true;
rank = false;
smoothness_factor = 0.5;
svGolay = true;
control_density = false;
matPatchColors = [0.75 0.75 0.75; 0.6 0.6 0.6; 0 0 1];
branch = zeros(1, numel(x));

for i=1:length(varargin)-1
    if(strcmp(varargin{i},'num_locs'))
        num_locs = varargin{i+1};
    elseif(strcmp(varargin{i},'avg_type'))
        avg_type = varargin{i+1};
    elseif(strcmp(varargin{i},'labels'))
        legend_flag = true;
        labels = varargin{i+1};
    elseif(strcmp(varargin{i},'show_error'))
        show_error = varargin{i+1};  
    elseif(strcmp(varargin{i},'normalize'))
        normalize = varargin{i+1};
    elseif(strcmp(varargin{i},'rank'))
        rank = varargin{i+1};
        if rank
            control_density = false;
        end
    elseif(strcmp(varargin{i},'svGolay'))
        svGolay = varargin{i+1};
    elseif(strcmp(varargin{i},'smooth'))
        smoothness_factor = varargin{i+1};
    elseif(strcmp(varargin{i},'branch'))
        branch = varargin{i+1};
    end
end

if (rank)
    x = tiedrank(x);
end    

weights = zeros(num_locs, length(x));
weights_win = zeros(num_locs, length(x));
Y_vals = zeros(num_locs, size(Y, 2));
Y_errs = zeros(num_locs, size(Y, 2));
Yn = size(Y, 1);

% check weights cache
hashMat = DataHash(x); 
hashVal = DataHash(string2hash([avg_type hashMat]) + smoothness_factor);

tic;
% Check for cached
[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
curr_path = [curr_path filesep];
cachefilename = [curr_path 'cachePlotAlongTimeResults.mat'];

%checking for old lnn results for the same data
fileCheck = exist (cachefilename, 'file');

% if no file
% if (fileCheck==0) 
%     % create new hash map
%     mapMat = containers.Map(); 
% else
%     try
%         % loading the old weights results
%         file = load(cachefilename); 
%         mapMat=file.mapMat;
%     catch
%         fileCheck=0;
%         % create new hash map
%         mapMat = containers.Map(); 
%     end
% end
% 
% % check for key in the hash map
% check = isKey(mapMat,hashVal); 
% 
% % if weights found in cache
% if (check==1) 
%     % no need to run lnn again->returning old result
%     value=values(mapMat,{hashVal});
%     W = value{1};
%     weights= W.weights;
%     weights_win = W.weights_win;
%     fprintf('weights loaded from cache: %gs\n', toc);
% else
    % compute a weight for each value (data point), at each plot location
    for i=1:num_locs
        weights(i, :) = compute_weights(x, (i/num_locs)*range(x)+min(x), avg_type, smoothness_factor);
%         weights_win(i, :) = compute_weights(x, (i/num_locs)*range(x)+min(x), 'sliding', smoothness_factor);
    end
    fprintf('weights computed: %gs\n', toc);

%     % while the hash map is too big removing the first element
%     while length(mapMat)>5 
%         lstKeys= keys(mapMat); 
%         remove(mapMat,lstKeys(1));
%     end
% 
%     % adding the name and lnn result to hashmap
%     W.weights = weights;
%     W.weights_win = weights_win;
%     mapMat(hashVal)= W; 
% end

tic;
%clean branch
if any(branch)
    pre_cleaning = weights;
    
    % locate trunk
    branchids = unique(branch)';
    maxidmean=0;
    trunkid =-1;
    for id=1:length(branchids)
        idmean=mean(weights(1:10,:)*(branch==branchids(id)));
        if idmean > maxidmean
            maxidmean = idmean;
            trunkid=id;
        end
    end
    
    % estimate branch point by quantity of branch flags along the trajectory 
    weights_by_branch = sum(weights(:, branch==trunkid), 2);
    weights_by_branch(:, end+1) = sum(weights(:, branch~=trunkid), 2);
    branching_loc = find(weights_by_branch(:, 1)<weights_by_branch(:, 2));
    branching_loc = branching_loc(1);
    
    % mute the effect of branch\trunk points beyond the branch point
    weights(branching_loc:num_locs, branch==trunkid) = 0;
    weights(1:branching_loc-1, branch~=trunkid) = 0;   
end
fprintf('estimating branch point: %gs\n', toc);


real_weights = weights;
for bri=1:2
    
    if any(branch)
        tic
        weights = real_weights;
        weights(:, branch==bri) = 0;
%         visible_locs = find(sum(weights, 2));
%         if (bri ~= trunkid) % if on a branch
%             trans_length = ceil(num_locs/20);
%             for transi = 1:(trans_length-1)
%                 weights(visible_locs(transi), branch~=bri) = (1-sqrt((transi/trans_length)))*pre_cleaning(visible_locs(transi), branch~=bri);
%             end
%         else
%             weights(visible_locs(end)+1, :) = pre_cleaning(visible_locs(end)+1, :);
%         end
        fprintf('correcting weights for transitioning: %gs\n', toc);
    end
    
    % Compute weighted averages at each location
    X = linspace(min(x), max(x), num_locs);
    Y_vals = weights*Y./repmat(sum(weights, 2), 1, size(Y, 2));

    if (svGolay)  
        for col=1:size(Y_vals, 2)
            Y_vals(:, col) = smooth(X, Y_vals(:, col),sqrt(num_locs*2), 'sgolay');
        end          
    end

    Y_vals_raw = Y_vals;
    if (normalize)
        % we want to normalize while accounting for branches
        y_vals_all = [];
        for ubri=unique(branch)'
            weights_tmp = real_weights;
            weights_tmp(:, branch~=ubri) = 0;
            % Compute weighted averages at each location
            y_vals_all = [y_vals_all; weights_tmp*Y./repmat(sum(weights_tmp, 2), 1, size(Y, 2))];
        end
        % we want to normalize to [0 1]
        mins = prctile(y_vals_all, 0, 1);
        Y_vals = bsxfun(@minus, Y_vals, mins);

        rngs = prctile(bsxfun(@minus, Y_vals_all, mins), 100, 1);
        Y_vals = bsxfun(@rdivide, Y_vals, rngs);

        Y_vals(Y_vals<0) = 0;        
        Y_vals(Y_vals>1) = 1;
    end

    matColors = distinguishable_colors(size(Y, 2));
    set(gca, 'ColorOrder', matColors);
    set(0, 'DefaultAxesColorOrder', matColors);

    marker = '-';
    if any(branch)
        if bri~=trunkid
            marker = Markers{bri};
        end
    end

    plot(X, Y_vals(:, 1),marker,...
         'LineWidth', 4,...
         'markersize', 6,...
         'Color', matColors(1, :)); 

    if (size(Y, 2)> 1)
        for col=2:size(Y, 2)
            hold on;
            plot(X, Y_vals(:, col),marker,...
             'LineWidth', 4,...
             'markersize', 6,...
             'Color', matColors(col, :));        
        end
    end
    
    if (show_error)

        % compute variace along X (symmetically)
        for i=1:num_locs       
            %symmetrical
            Y_errs = bsxfun(@minus,Y,(Y_vals_raw(i, :)));

            M = sum(weights(i, :)~=0);
            s = (M-1)/M;
            w_sum = sum(weights(i, :));

            Y_valerrs(i, :) = sqrt((weights(i, :)*((Y_errs).^2))/(s*w_sum));
        end

        if (normalize)
            Y_valerrs = bsxfun(@rdivide,Y_valerrs,rngs);
        end

        % plot the variance as a pretty translucent cloud around line
        for yi=1:size(Y, 2)

            % plot light grey background first for variance
            fill( [X, fliplr(X)], [(Y_vals(:, yi)-Y_valerrs(:, yi)./2)', fliplr((Y_vals(:, yi)+Y_valerrs(:, yi)./2)')],...
                matColors(yi,:),'linestyle','none','facealpha',.5);
        end
    end

    % Hold on if plotting more branch lines
    if (any(branch))
        hold on;
    end

end

    % show density histogram under the plot to show the concentration 
    if false %~(rank || control_density)    
        try
            dens = sum(weights_win, 2)';
            ca = axis;
            hold on;
            y_range = ca(4)-ca(3);
            y_buffer = y_range*.04;
        %     ylim([ca(3)-.1 ca(4)]);
            imagesc(X, [ca(3)-y_buffer ca(3)-(.5*y_buffer)], dens, [0, max(dens)]);   
            colorbar;
            axis([ca(1:2) ca(3)-y_buffer ca(4)]);
        catch e
            disp(getReport(e,'extended'));
        end
    end

if (legend_flag)
    legend(remove_repeating_strings(labels), 'Interpreter', 'none');
end

% if (check==0)
%     % saving into file
%     try
%     save(cachefilename,'mapMat');
%     catch
%         fprintf('error caching weights in %s', cachefilename);
%     end
% end

end

function weights = compute_weights(points, loc, type, factor)
    
    range = quantile(points, .98) - quantile(points, .02);
    min_std_dev = factor*.1*range; % minimum std_dev for dense regions
    max_std_dev = .19*range; % max std_dev for sparse regions  
    linear_slope = 10/range;
    
    if strcmpi('sliding', type) %set '1's on the indices in the windows 
        weights = (points < (loc + 2*(min_std_dev))) & ...
                  (points > (loc - 2*(min_std_dev)));
    
    elseif strcmpi('linear', type)
        weights = 1 - linear_slope*(abs(points - loc));
        weights(weights<0) = 0;
        weights = weights/sum(weights);       
    
    elseif strcmpi('squared', type)
        weights = 1 - ((linear_slope*(points - loc)).^2);
        weights(weights<0) = 0;
        weights = weights/sum(weights);
             
    elseif strcmpi('gaussian_var', type)
        [f, xi] = persistent_ksdensity(points);
    	d = f(minind(abs(xi-loc)))/max(f);

        std_dev = d*(min_std_dev)+(1-d)*max_std_dev;

        weights = ((2*pi*(std_dev^2))^(-1))*exp(-.5*((points - loc)/std_dev).^2);
    
    else % default is strcmpi('gaussian', type)
        weights = ((2*pi*(min_std_dev)^2)^(-1))*exp(-.5*((points - loc)/min_std_dev).^2);
    end
end

function [f, xi] = persistent_ksdensity(points)

    persistent f_p;
    persistent xi_p;    

    if isempty(f_p)
        [f_p, xi_p] = ksdensity(points);
    end
    f = f_p;
    xi = xi_p;
end