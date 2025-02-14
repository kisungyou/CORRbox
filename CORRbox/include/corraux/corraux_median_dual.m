function [logXY, dval] = corraux_median_dual(X,Y)
    % compute log and distance at the same time
    D = corraux_findD(X, Y);
    Z = D*Y*D;
    logXY = corraux_log(X, Z);
    dval  = corraux_spddist(X, Z);
end