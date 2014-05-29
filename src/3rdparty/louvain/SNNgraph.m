function G = SNNgraph( x, k, distance )
% G = SNNgraph( x, k )
% Local version. Analogous server version: large_SNN_graph.m
% Compute k-nearest neighbors for each point in x, broken into 1000-point
% chunks for efficiency, then compute a sparse secondary-nearest-neighbor 
% (SNN) graph with simcos as weights (#{intersection}/k)

n = size( x, 1 );
if n > 6000
    chunk_size = 1000;
else
    chunk_size = n;
end

if nargin < 2
    k = 30;
end

if nargin < 3
	distance = 'euclidean'
end

iters = ceil( n / chunk_size );
iter = 1;

IDX = NaN( n, k );
fprintf(1, 'Computing %i nearest neighbors\n', k)
for from = 1:chunk_size:n
    t = tic;
    
    to = min( from + chunk_size - 1, n );
    rx = (from:to)';
    x_from_to = x( rx, : );
    
    idx = knnsearch( x, x_from_to, 'k', k+1, 'distance', distance );
    IDX(rx,:) = idx(:,2:end);
    
    fprintf(1, 'Iteration %i of %i complete: %.1f seconds\n', iter, iters, toc(t) )
    iter = iter + 1;
end
    
fprintf(1, 'Building SNN graph...\n' )
I = NaN( 1, k*n );
J = I;
S = I;
row = 1;
pctdone = 10;
for i = 1:n
    
    neighbors = IDX(i,:);
%     neighbors( neighbors == i ) = [];
    neighbs_of_neighbs = IDX(neighbors,:);
%     neighbs_of_neighbs(:,1) = [];
    intersection_size = arrayfun(@(x) length(intersect(neighbors, neighbs_of_neighbs(x,:))), ...
        1:k);
    simcos = intersection_size ./ k;
%     I = [I; repmat(i,1,k)];
%     J = [J; neighbors];
%     S = [S; simcos];
    idx = row:row+k-1;
    I(idx) = repmat(i,1,k);
    J(idx) = neighbors;
    S(idx) = simcos;
    row = row + k;
   
    if ~mod(i,n/10)
        fprintf(1, '%i percent complete\n',pctdone)
        pctdone = pctdone + 10;
    end
end

G = sparse( I, J, S, n, n );