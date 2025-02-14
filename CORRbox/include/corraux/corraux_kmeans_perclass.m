function [centers, variation] = corraux_kmeans_perclass(input, label, geometry)

% given a label vector, compute K cluster centers with cluster-specific
% variation so that this information can be used in the main kmeans
% algorithm. 

% compute mean per class
ulabel = sort(unique(label));
K = length(ulabel);
p = size(input.data, 1);
geomobj = corrbox_geomhandle(geometry);

% empty arrays for recording
centers = zeros(p,p,K);
variation = zeros(K,1);

% iterate per class
for i=1:K
    idx = (label==ulabel(i));
    [part_mat, part_var] = geomobj.fmean(input.data(:,:,idx));
    centers(:,:,i) = part_mat;
    variation(i)   = part_var;
end
end

