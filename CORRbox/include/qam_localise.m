function output = qam_localise(Cloc, C3d)

% localize the data near the given point.
% QAM : take logarithmic map to the location, and return the tril.

% parameteres
[p,~,N] = size(C3d);

% prepare a return object
output = zeros(N,round(p*(p-1)/2));


corraux_log(frechet_mean, C3d(:,:,n));

% iterate
for n=1:N
    tgt = corraux_log(Cloc, C3d(:,:,n));
    output(n,:) = trf_tril_to_vec(tgt, false);
end
end