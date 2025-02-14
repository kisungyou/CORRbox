function LT1mat = trf_COR_to_LT1(Cmat)

% input
%   Cmat   : correlation matrix
% output
%   LT1mat : lower-triangular with unit diagonals

%% main 
L = chol(Cmat, "lower");
LT1mat = diag(1./diag(L))*L;

end