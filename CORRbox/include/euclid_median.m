function [Cmed, Cvar] = euclid_median(C3d, maxiter, thr)

% naive euclidean geometry + frechet median. working as a baseline. 
% for this, implement Euclidean weiszfeld algorithm as of corraux variables. 

% case 1 : just a single observation
if length(size(C3d)) < 3
    Cmed = C3d;
    Cvar = 0;
% case 2 : multiple observations
else
    % Frechet median
    Cmed = corraux_weiszfeld(C3d, maxiter, thr);

    % compute the variation
    N = size(C3d, 3);
    varsq = zeros(1,N);
    for n=1:N
        varsq(n) = euclid_dist(Cmed, C3d(:,:,n));
    end
    Cvar = mean(varsq);
end
end