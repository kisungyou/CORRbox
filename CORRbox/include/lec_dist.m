function distval = lec_dist(C1, C2)

% given two correlation matrices, compute the distance between two
% correlation matrices under the LEC (Log-Euclidean-Cholesky) geometry.

%% main
%  conversion
LogL1 = trf_COR_to_LogLT1(C1);
LogL2 = trf_COR_to_LogLT1(C2);

%  compute the norm
distval = norm(LogL1-LogL2, "fro");

end