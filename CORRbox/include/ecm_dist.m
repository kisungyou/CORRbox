function distval = ecm_dist(C1, C2)

% given two correlation matrices, compute the distance between two
% correlation matrices under the ECM (Euclidean-Cholesky Metric) geometry.

%% main
%  convert into LT1 vectorial form
L1vec   = aux_chol2vec(trf_COR_to_LT1(C1));
L2vec   = aux_chol2vec(trf_COR_to_LT1(C2));

%  simple Euclidean distance
distval = norm(L1vec-L2vec, 2);

end