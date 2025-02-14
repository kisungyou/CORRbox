function [embedding, Cmean, apply_mean, apply_proj] = common_pga_main(C3d, ndim, geometry)

% find the coordinates around the mean
gobj  = corrbox_geomhandle(geometry);
[Cmean,~] = gobj.fmean(C3d);
X_train = gobj.localise(C3d, Cmean);
[N,P] = size(X_train);

% compute the mean and normalize
X_scale = zeros(N,P);
apply_mean = mean(X_train, 1);
for n=1:N
    X_scale(n,:) = X_train(n,:) - apply_mean;
end

% economical SVD
[~,~,V] = svd(X_scale, "econ");
apply_proj = V(:,1:round(ndim));

% embedding
embedding = X_scale*apply_proj;

end