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

Markers={'-', '--', ':'};

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
branchY = zeros(1, numel(x));

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
    elseif(strcmp(varargin{i},'branchY'))
        branchY = varargin{i+1};
        Y_scale = branchY-median(branchY);
        Y_scale(Y_scale<0) = Y_scale(Y_scale<0)./max(abs((Y_scale(Y_scale<0))));
        Y_scale(Y_scale>0) = Y_scale(Y_scale>0)./max(Y_scale(Y_scale>0));

    end
end

if (rank)
    x = tiedrank(x);
end    

real_weights = zeros(num_locs, length(x));
weights_win = zeros(num_locs, length(x));
Y_vals = zeros(num_locs, size(Y, 2));
Y_errs = zeros(num_locs, size(Y, 2));
Yn = size(Y, 1);

tic;

% compute a weight for each value (data point), at each plot location
for i=1:num_locs
    real_weights(i, :) = compute_weights(x, (i/num_locs)*range(x)+min(x), avg_type, smoothness_factor);
end
fprintf('weights computed: %gs\n', toc);

tic

% Compute weighted averages at each location
y_vals_all = real_weights*Y./repmat(sum(real_weights, 2), 1, size(Y, 2));
fprintf('Eighted marker values computed: %gs\n', toc);

tic
if normalize
    % we want to normalize to [0 1]
    mins = prctile(y_vals_all, 0, 1);
    rngs = prctile(bsxfun(@minus, y_vals_all, mins), 100, 1); 
    
    y_vals_all = bsxfun(@minus, y_vals_all, mins);
    y_vals_all = bsxfun(@rdivide, y_vals_all, rngs);

    y_vals_all(y_vals_all<0) = 0;        
    y_vals_all(y_vals_all>1) = 1;
end
fprintf('Normalization values computed: %gs\n', toc);
    
Y_vals_branches = cell(1,2);
Y_vals_raws     = cell(1,2);

% compute both sides of wine glass
for bri=1:2
    Y_scale = -1 * Y_scale;
    
    tic
    weights = real_weights;
    weights=weights.*repmat(1-(max(Y_scale(:)', 0)).^.005, num_locs, 1);
    fprintf('correcting weights for transitioning: %gs\n', toc);
    
    % Compute weighted averages at each location
    X = linspace(min(x), max(x), num_locs);
    
    Y_vals_raws{bri} = weights*Y./repmat(sum(weights, 2), 1, size(Y, 2));    

    if (normalize)       
        Y_vals = bsxfun(@minus, Y_vals_raws{bri}, mins);
        Y_vals = bsxfun(@rdivide, Y_vals, rngs);

        Y_vals(Y_vals<0) = 0;        
        Y_vals(Y_vals>1) = 1;
    end
    
    Y_vals_branches{bri} = Y_vals;    
end

Y_vals_main = Y_vals_branches{1};
Y_vals = Y_vals_branches{2};

% branch line - we are selective on the plotting
for loc=num_locs:-1:2
    markers = find(abs((Y_vals(loc-10,:) - Y_vals_main(loc-10, :))) < .25 & abs((Y_vals(loc,:) - Y_vals_main(loc, :))) < .3);
    Y_vals(1:(loc-16),markers) = nan;
%     Y_vals_main(1:(loc-20),markers) = y_vals_all(1:(loc-20),markers);

    if all(isnan(Y_vals(loc,:)))
        markers = isnan(Y_vals(num_locs-21,:));
        Y_vals(:,markers) = nan;
        break;
    end        
end
Y_vals_main(isnan(Y_vals)) = y_vals_all(isnan(Y_vals))


if (svGolay)  
    for col=1:size(Y_vals, 2)
        Y_vals_main(:, col) = smooth(X, Y_vals_main(:, col),sqrt(num_locs*2), 'sgolay');
    end          
end

Y_vals_branches{1} = Y_vals_main;
Y_vals_branches{2} = Y_vals;

% iterate for plotting
for bri=1:2
    matColors = distinguishable_colors(size(Y, 2));
    set(gca, 'ColorOrder', matColors);
    set(0, 'DefaultAxesColorOrder', matColors);

    % change marker selection
    marker = Markers{bri};

    Y_vals = Y_vals_branches{bri};
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
        Y_vals_raw = Y_vals_raws{bri};

        % compute variace along X (symmetically)
        for i=1:num_locs       
            
            %symmetrical
            Y_errs = bsxfun(@minus,Y,(Y_vals_raw(i, :)));

            M = sum(weights(i, :)~=0);
            s = (M-1)/M;
            w_sum = sum(weights(i, :));

            Y_valerrs(i, :) = .5*sqrt((weights(i, :)*((Y_errs).^2))/(s*w_sum));
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
	hold on;
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
    l=legend(remove_repeating_strings(labels), 'Interpreter', 'none');
    set(l, 'Location','NorthEastOutside');
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