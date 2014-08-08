function clusterHeatmap( data, labels, names, normalize )
% clusterHeatmap( data, labels, names, normalize )
% Simplified code; old version is commented out
% Any string entered as the forth argument selects the option to scale each
% column of data between 0 and 1
% Names can be skipped by supplying an empty matrix

if size(data,1) ~= length(labels)
    error('data and labels are not the same length')
end

if nargin > 2
    if ~isempty(names)
        namesgiven = true;
    end
else
    namesgiven = false;
end

if nargin > 3
    normalize = true;
else
    normalize = false;
end

data( labels == 0, :) = [];
labels( labels == 0 ) = [];
clusterSize = arrayfun(@(x) sum(labels==x), 1:max(labels) );

if normalize
    % scale each column to [0 1]
    % apply ceiling at 99th percentile
    mtx = NaN( size(data) );
    for i = 1:size(data,2)
        cutoff = prctile(data(:,i), 99);
        v = data(:,i) ./ cutoff;
        v( v > 1 ) = 1;
        mtx(:,i) = v;
    end
    
else
    
    mtx = data;
    
end

[~,ix] = sort(labels);
mtx = mtx(ix,:)';
lines = cumsum( clusterSize ) + .5;
%imagesc(mtx, [-max(max(abs(mtx))),max(max(abs(mtx)))]);    %use this
%option with red/blue heatmap
imagesc(mtx);
colormap(linspecer);
%colormap(interpolate_colormap(redbluecmap,64));    %use this for red/blue
%heatmap centered around zero
for i = 1:length(lines)-1 % don't plot first or last line
    hold on;
    plot( [lines(i) lines(i)], ylim(), 'k--', 'LineWidth', 2);
    hold off;
end

if namesgiven
    set(gca,'ytick',1:length(names))
    set(gca,'yticklabel', names, 'fontsize', 14)
end

% remove empty clusters
remove = clusterSize==0;
lines(remove) = [];
k = length(unique(labels));
% set xticks in the middle of each cluster
midpoints = [];
lines = [0 lines];
for i = 1:k
    midpoints(i) = (lines(i) + lines(i+1))/2;
end
set(gca,'xtick',midpoints, 'xticklabel', unique(labels))
xlabel('cluster ID', 'fontsize', 14 )

% 
% function argOut = clusterHeatmap( data, labels, varargin )
% % argOut = clusterHeatmap( data, labels, varargin )
% % OPTIONS
% % 'sortrows', t/f: Sort rows by hierarchical clustering of cluster means
% %   default: false
% % 'plotmedians', t/f: Instead of individual observations, plot cluster
% %   medians
% %   default:false
% % 'names', cellarray of strings: names for each column of data (row of
% %   heatmap)
% sortrows = false;
% plotmedians = false;
% channelnames = {};
% means = [];
% medians = [];
% roworder = [];
% samplePercent = [];
% 
% for i = 1:2:length(varargin)
%     
%     switch lower(varargin{i})
%         
%         case 'sortrows'
%             
%             sortrows = varargin{i+1};
%             
%         case 'plotmedians'
%             
%             plotmedians = varargin{i+1};
%             
%         case 'channelnames'
%             
%             channelnames = varargin{i+1};
%             
%         case 'names' 
%             
%             channelnames = varargin{i+1};
%             
%     end
% end
%             
% data( labels == 0, :) = [];
% labels( labels == 0 ) = [];
% clusterSize = arrayfun(@(x) sum(labels==x), 1:max(labels) );
% 
% % sort clusters by similarity
% if sortrows
%     means = [];
%     for c = 1:max(labels)
%         means( c, : ) = mean( data( labels == c, : ) );
%     end
%     means = means';
%     rowdist = pdist( means, 'seuclidean' );
%     rowtree = linkage( rowdist, 'average' );
%     roworder = optimalleaforder( rowtree, rowdist );
% else
%     roworder = 1:size(data,2);
% end
% 
% % figure
% 
% if plotmedians
%     
%     samplePercent = 100*clusterSize ./ sum(clusterSize);
%     medians = [];
%     for c = 1:max(labels)
%         medians( c, : ) = median( data( labels == c, : ) );
%     end
%     medians = medians';
%     if isempty(channelnames)
%         error('Plotmedians requires rownames')
%     end
%     heatmap( medians( roworder, : ), samplePercent, channelnames(roworder), [], ...
%         'GridLines', '-', 'ShowAllTicks', true, 'TickAngle', 45 );
%     
% else
%     
%     [~,ix] = sort(labels);
%     mtx = data';
%     mtx = mtx( roworder, ix );
%     lines = cumsum( clusterSize ) + .5;
%     imagesc(mtx);
%     for i = 1:length(lines)
%         hold on;
%         plot( [lines(i) lines(i)], ylim(), 'k--', 'LineWidth', 2 );
%         hold off;
%     end
%     if ~isempty( channelnames )
%         set(gca, 'ytick', 1:length(channelnames))
%         set(gca, 'yticklabel', channelnames(roworder) )
%     end
%     
% end
% 
% if ~isempty( medians )
%     argOut.medians = medians( roworder, : );
% %     argOut.means = means( roworder, : );
%     argOut.percents = samplePercent;
%     argOut.roworder = roworder;
%     argOut.rowNames = channelnames(roworder);
% else
%     argOut.percents = samplePercent;
% end