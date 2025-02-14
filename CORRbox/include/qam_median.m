function [Cmed, Cvar] = qam_median(C3d, maxiter, thr)

% revisiting quotient affine geometry of correlation manifold for computing
% the Frechet median.

% case 1 : just a single observation
if length(size(C3d)) < 3
    Cmed = C3d;
    Cvar = 0;
    % case 2 : multiple observations
else
    % step 1. initialize
    p = size(C3d, 1);
    N = size(C3d, 3);
    median_old = mean(C3d, 3); % ambient average

    % step 2. iterate
    vec_dist = zeros(1,N);
    vec_log  = zeros(p,p,N);

    for it=1:maxiter
        % 1. compute log and distance
        for n=1:N
            [tmpLogmat, tmpdist] = corraux_median_dual(median_old, C3d(:,:,n));
            vec_log(:,:,n) = tmpLogmat;
            vec_dist(n)    = tmpdist;
        end

        % 2. updating information
        tmpLog = zeros(p,p);
        tmpNum = 0;
        for n=1:N
            if (vec_dist(n) > 1e-10)
                tmpNum = tmpNum + (1/vec_dist(n));
                tmpLog = tmpLog + (vec_log(:,:,n)/vec_dist(n));
            end
        end

        % 3. do the update
        median_new = corraux_exp(median_old, tmpLog/tmpNum, 1.0);
        increment  = norm(median_old-median_new,"fro");
        median_old = median_new;

        % 4. stopping criterion & update
        if (increment < thr)
            break;
        end
    end

    % compute Frechet variation
    variation_vec = zeros(1,N);
    for i=1:N
        variation_vec(i) = corraux_dist(median_old, C3d(:,:,i));
    end
    variation = mean(variation_vec);

    % wrap it
    Cmed = median_old;
    Cvar = variation;
end
end

