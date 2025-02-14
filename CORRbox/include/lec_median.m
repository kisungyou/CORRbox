function [Cmed, Cvar] = lec_median(C3d, maxiter, thr)

% lec geometry

% case 1 : just a single observation
if length(size(C3d)) < 3
    Cmed = C3d;
    Cvar = 0;
% case 2 : multiple observations
else
    % Extrinsic 
    p = size(C3d, 1);
    N = size(C3d, 3);
    trfC3d = zeros(p,p,N);
    for n=1:N
        trfC3d(:,:,n) = trf_COR_to_LogLT1(C3d(:,:,n));
    end

    % Frechet median
    trfCmed = corraux_weiszfeld(trfC3d, maxiter, thr);
    Cmed = trf_LogLT1_to_COR(trfCmed);
    Cmed = (Cmed + Cmed')/2;

    % compute the variation
    N = size(C3d, 3);
    varsq = zeros(1,N);
    for n=1:N
        varsq(n) = lec_dist(Cmed, C3d(:,:,n));
    end
    Cvar = mean(varsq);
end
end