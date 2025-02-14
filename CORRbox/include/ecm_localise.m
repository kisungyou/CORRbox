function output = ecm_localise(Cloc, C3d)

% localize the data near the given point.
% ECM : apply transformation, subtract the location, and return the tril.

% parameteres
[p,~,N] = size(C3d);

% prepare a return object
output = zeros(N,round(p*(p-1)/2));
trfloc = trf_COR_to_LT1(Cloc);

% iterate
for n=1:N
    tgt = trf_COR_to_LT1(C3d(:,:,n))-trfloc;
    output(n,:) = trf_tril_to_vec(tgt, false);
end
end