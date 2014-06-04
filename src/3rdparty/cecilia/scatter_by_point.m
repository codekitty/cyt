function scatter_by_point(x, y, color_var, dot_size)
    
    %figure;
    x_lim = [min(x)-min(x)/20,max(x)+max(x)/20];
    y_lim = [min(y)-min(y)/20, max(y)+max(y)/20];
    
    x_lim = [min(x)-(max(x)-min(x))*0.03, max(x)+(max(x)-min(x))*0.03];
    y_lim = [min(y)-(max(y)-min(y))*0.03, max(y)+(max(y)-min(y))*0.03];
    
    colors = distinguishable_colors(max(color_var));
    
    for i=1:size(x),

        scatter(x(i),y(i), dot_size(i), colors(color_var(i),:), 'fill');

        xlim(x_lim);
        ylim(y_lim);
        hold all;
        
    end
    
    %xlabel('tSNE1');
    %ylabel('tSNE2');
    
    
    %adding legend
    legend_handle = legend(gca,num2str(unique(color_var)), 'location', 'eastoutside');
    legend_markers = findobj(get(legend_handle, 'Children'), 'marker', 'o');    %find handles of markers in legend
    
    p = length(legend_markers);
    for i=1:length(legend_markers), %set color of markers in legend
        set(legend_markers(p), 'markerfacecolor', colors(i,:));
        p=p-1;
    end

end