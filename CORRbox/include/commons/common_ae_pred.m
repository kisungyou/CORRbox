function embedding = common_ae_pred(newcorr, geometry, Cmean, aemodel)

% just check an input
if (~corraux_checker(newcorr))
    error("* corr_ae : the new input is not a valid object. Use 'corr_init'.");
end

% localize the new data
gobj   = corrbox_geomhandle(geometry);
X_test = gobj.localise(newcorr.data, Cmean);
X_test = X_test';

% apply the encoder
embedding = encode(aemodel, X_test);
embedding = embedding'; % convert to row majors

end