function [embedding, pred] = common_ae_main(C3d, ndim, geometry, alpha, beta)

% find the coordinates around the mean
gobj  = corrbox_geomhandle(geometry);
[Cmean,~] = gobj.fmean(C3d);
X_train = gobj.localise(C3d, Cmean);
X_train = X_train'; % trainAutoencoder takes column-major data

% train the autoencoder
aemodel = trainAutoencoder(X_train, ndim, ...
    "L2WeightRegularization",alpha, ...
    "SparsityRegularization",beta, ...
    "ShowProgressWindow", false, ...
    "ScaleData", false);

% compute the embedding using the encoder
embedding = encode(aemodel, X_train);
embedding = embedding'; % change to row-major

% predict function
pred = @(newobj) common_ae_pred(newobj, geometry, Cmean, aemodel);

end
