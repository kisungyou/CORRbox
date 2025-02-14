function embedding = corraux_PCA(X, ndim)

% this is my naive implementation of PCA using the standard idea. we assume
% that the matrix is following the convention of (N x P) format.

% data normalization
[N,P] = size(X);
meanX = mean(X,1);
Y = zeros(N,P);
for n=1:N
    Y(n,:) = X(n,:)-meanX;
end

% economical SVD
[~,~,V] = svd(Y, "econ");
proj = V(:,1:round(ndim));

% projection
embedding = Y*proj;

end