function embedding = ecm_pga(C3d, ndim)

% PGA with ECM metric

% find the coordinates around the mean
highD = ecm_coord_near_mean(C3d);
if (size(highD,2) >= ndim)
    error("* corr_pga : the provided 'ndim' should be smaller.");
end
% apply PCA
embedding = corraux_PCA(highD, round(ndim));

end