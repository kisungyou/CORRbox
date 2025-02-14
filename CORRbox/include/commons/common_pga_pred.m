function embed_new = common_pga_pred(newobj, Cmean, apply_mean, apply_proj, geometry)

% just check an input
if (~corraux_checker(newobj))
    error("* corr_pga : the new input is not a valid object. Use 'corr_init'.");
end

% localise again
gobj   = corrbox_geomhandle(geometry);
X_test = gobj.localise(newobj.data, Cmean);
[N,P] = size(X_test);
if (length(apply_mean)~=P)
    error("* corr_pga : the new input does not match the dimensionality of the trained model.");
end

% normalize + embedding
X_scale = zeros(N,P);
for n=1:N
    X_scale(n,:) = X_test(n,:)-apply_mean;
end
embed_new = X_scale*apply_proj;

end