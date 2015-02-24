function [centroids, cluster_mapping,cluster_sizes,cellsInCluster] = compute_cluster_centroids(data, inds, cluster_channel)


    cluster_mapping = [];
    centroids = zeros(0,size(data, 2));
    k = 1;
    cluster_sizes = [];
    
    for i=1:size(inds,2), %looping through selected gates (samples)
        cluster_mapping = vertcat(cluster_mapping, horzcat(repmat(i,[size(unique(cluster_channel(inds{i}))),1]),unique(cluster_channel(inds{i}))));
        unique_clusters = unique(cluster_channel(inds{i}));
        
        for j=1:length(unique_clusters), %looping through clusters in gate/sample
            sub_data = data(inds{i},:); %getting the cells of the choosen gates
            centroids(end+1,:) = mean(sub_data(cluster_channel(inds{i}) == unique_clusters(j),:));
            cellsInCluster(k) = size(data(cluster_channel(inds{i}) == unique_clusters(j),:),1);
            cluster_sizes(k) = cellsInCluster(k)/ size(cluster_channel(inds{i}),1);
            %cluster_sizes(k) = size(data(cluster_channel(inds{i}) == unique_clusters(j),:),1) / size(cluster_channel(inds{i}),1);
            k=k+1;
        end
    end
    
    cluster_mapping = horzcat(linspace(1,size(cluster_mapping,1),size(cluster_mapping,1))', cluster_mapping, cluster_sizes');   %col1 = cluster number, col2=sample number, col3=withing sample cluster number
    
    
    % mapping between clusters and meta clusters
    meta_clusters = [];
    
    for i=1:size(inds,2), %looping through samples
        unique_clusters = unique(cluster_channel(inds{i}));
        meta = cluster_channel(inds{i});
        
        for j=1:length(unique_clusters),  %looping through unique clusters in sample
            meta_clusters = vertcat(meta_clusters, meta(find(cluster_channel(inds{i}) == unique_clusters(j), 1,'first')));
        end
    end
    
    cluster_mapping = horzcat(meta_clusters, cluster_mapping);

    %cleaning the output
    cluster_mapping(isnan(centroids(:,1)),:) = [];    %removing rows that correspond to clusters where centroid is NaN
    centroids(isnan(centroids(:,1)),:) = [];  %removing rows with NaNs in centroids
    centroids=centroids(find(cluster_mapping(:,1)),:); %removing rows with 0 in cluster number
    cluster_mapping=cluster_mapping(find(cluster_mapping(:,1)),:);
    
    
end