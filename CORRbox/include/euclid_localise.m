function output = euclid_localise(Cloc, C3d)

% localize the data near the given point.
% EUCLIDEAN : subtract the location and return the lower triangular part.

% parameteres
[p,~,N] = size(C3d);

% prepare a return object
output = zeros(N,round(p*(p-1)/2));

% iterate
for n=1:N
    tgt = C3d(:,:,n)-Cloc;
    output(n,:) = trf_tril_to_vec(tgt, false);
end
end