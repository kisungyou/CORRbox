function [Cmean, Cvar] = euclid_mean(C3d)

% given a 3d array (p,p,N) of correlation matrices, compute the frechet
% mean and variation under the Euclidean geometry.
% if a single matrix is given, just return a matrix with variation 0

if length(size(C3d)) < 3
    Cmean = C3d;
    Cvar = 0;
else
    %% main
    %  get the size of an array
    N = size(C3d, 3);

    %  compute the mean
    Cmean = mean(C3d, 3);

    %  compute the Frechet variation
    Cvar = 0;
    for n=1:N
        Cvar = Cvar + (1/N)*(norm(Cmean-C3d(:,:,n), "fro")^2);
    end
end
end