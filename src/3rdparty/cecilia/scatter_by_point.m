function scatter_by_point(x, y, gcolors, dot_size)
        
    groups = unique(gcolors);
    nGroups = numel(groups);
    
    %color_per_group = distinguishable_colors(nGroups);
    
    %only 24 colors
    map=[[0,0.447,0.741];    [0.85,0.325,0.098]; [0.929,0.694,0.125];
         [0.494,0.184,0.556];[0.466,0.674,0.188];[0.201,0.845,0.933];
         [0,1,1];            [1,0,0];            [0,1,0];
         [0.194,0.684,0.256];[0.166,0.374,0.088];[0.101,0.545,0.733];
         [0,0.247,0.541];    [0.65,0.125,0.098]; [0.729,0.494,0.025];
         [0,0,1];            [1,1,0];            [1,0,1]
         [0.2,0.647,0.941];  [0.95,0.525,0.298]; [0.929,0.894,0.525];
         [0.635,0.078,0.184];[0.766,0.974,0.388];[0.201,0.945,0.533]];
     
    color_per_group = colormap(map);
    
    % plot each color in a seperate group (otherwise the legend won't work)
    for i=1:nGroups
        curr_group = groups(i);
        curr_color_inds = (curr_group==gcolors);

        scatter(x(curr_color_inds),y(curr_color_inds),...
                dot_size(curr_color_inds),...
                color_per_group(i, :), 'fill');

        hold on;      
    end
    
	% find axis limits
    x_range  = max(x)-min(x);
    x_buffer = x_range*.03;
    x_lim    = [min(x)-x_buffer, max(x)+x_buffer];
    
    y_range  = max(y)-min(y);
    y_buffer = y_range*.03;
    y_lim    = [min(y)-y_buffer, max(y)+y_buffer];

    xlim(x_lim);
    ylim(y_lim);

end