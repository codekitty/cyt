function scatter_by_point(x, y, gcolors, dot_size)
        
    groups = unique(gcolors);
    nGroups = numel(groups);
    
    color_per_group = distinguishable_colors(nGroups);
    
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