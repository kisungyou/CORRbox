function [Cmean, Cvar] = qam_mean(C3d, maxiter, thr)

% given a 3d array (p,p,N) of correlation matrices, compute the frechet
% mean and variation under the Quotient Affine Geometry.
% if a single matrix is given, just return a matrix with variation 0

% control stopping criteria
if nargin < 2
    maxiter = 100;
    thr = 1e-8;
end

if length(size(C3d)) < 3
    Cmean = C3d;
    Cvar  = 0;
else
    %% MAIN COMPUTATION
    %  preprocessing
    [p,~,N] = size(C3d);

    %  initialization in the ambient space
    mean_old = mean(C3d, 3);

    %  iterate
    increment = 10000;
    tvecs     = zeros(p,p,N);
    iter      = 1;
    while (increment > thr)
        %   2-1. compute logarithmic map
        for i=1:N
            tvecs(:,:,i) = corraux_log(mean_old, C3d(:,:,i));
        end
        %   2-2. mean
        mean_tvecs = mean(tvecs,3);
        %   2-3. update using exponential map
        mean_new   = corraux_exp(mean_old, mean_tvecs);
        %   2-4. compute increment
        increment  = corraux_dist(mean_old, mean_new);
        %   2-5. update old as new
        mean_old   = mean_new;
        %   2-6. iteration count
        iter = iter + 1;
        if (iter >= maxiter)
            break;
        end
    end

    %  step 3. compute Frechet variation
    Cmean         = mean_old;
    variation_vec = zeros(1,N);
    for i=1:N
        variation_vec(i) = corraux_dist(Cmean, C3d(:,:,i))^2;
    end
    Cvar = mean(variation_vec);
end
end