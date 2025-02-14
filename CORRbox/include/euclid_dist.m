function distval = euclid_dist(C1, C2)

% given two correlation matrices, compute the distance between two
% correlation matrices under the ambient geometry of Euclidean space.

distval = norm(C1-C2, "fro");

end