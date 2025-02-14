function LT1log = trf_COR_to_LogLT1(Cmat)

% input
%   Cmat   : correlation matrix
% output
%   LT1log : lower-triangular with unit diagonals with logm applied

%% main 
L = chol(Cmat, "lower");
LT1mat = diag(1./diag(L))*L;
LT1log = logm(LT1mat);

end