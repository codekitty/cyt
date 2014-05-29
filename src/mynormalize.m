function data = mynormalize(data, percentile)
% normalize data according to specified percentile

            disp('Normalizing according to the 99th percentile...');
            data = data-repmat(prctile(data, 100-percentile, 1), size(data,1),1);
            data = data./repmat(prctile((data), percentile, 1),size(data,1),1);

            data(data > 1) = 1;
            data(data < 0) = 0;
            data(isinf(data)) = 0;
end