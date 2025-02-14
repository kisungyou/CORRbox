function distmat = euclid_pdist(C3d)

% size
N = size(C3d, 3);

% prepare
distmat = zeros(N,N);

% compute
for i=1:(N-1)
    for j=(i+1):N
        dval = euclid_dist(C3d(:,:,i), C3d(:,:,j));

        distmat(i,j) = dval;
        distmat(j,i) = dval;
    end
end
end