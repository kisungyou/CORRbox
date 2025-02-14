function [Cmean, Cvar] = lec_mean(C3d)

% given a 3d array (p,p,N) of correlation matrices, compute the frechet
% mean and variation under the LEC geometry.
% if a single matrix is given, just return a matrix with variation 0

if length(size(C3d)) < 3
    Cmean = C3d;
    Cvar  = 0;
else
    %% main
    %  get the size of an array
    [p,~,N] = size(C3d);

    %  create an empty array
    Ctmp = zeros(p,p);

    %  convert each object and sum it
    for n=1:N
        Ctmp = Ctmp + (1/N)*trf_COR_to_LogLT1(C3d(:,:,n));
    end

    %  inverse transform + symmetrization
    Cmean = trf_LogLT1_to_COR(Ctmp);
    Cmean = (Cmean + Cmean')/2;

    %  compute Frechet variation
    Cvar = 0;
    for n=1:N
        Cvar = Cvar + (1/N)*(lec_dist(Cmean, C3d(:,:,n))^2);
    end
end
end