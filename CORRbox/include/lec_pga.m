function embedding = lec_pga(C3d, ndim)

% PGA with LEC metric

% find the coordinates around the mean
highD = lec_coord_near_mean(C3d);
if (size(highD,2) >= ndim)
    error("* corr_pga : the provided 'ndim' should be smaller.");
end
% apply PCA
embedding = corraux_PCA(highD, round(ndim));

end