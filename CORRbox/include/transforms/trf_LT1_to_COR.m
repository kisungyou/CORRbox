function Cmat = trf_LT1_to_COR(LT1mat)

% input
%   LT1mat : lower-triangular with unit diagonals
% output
%   Cmat   : correlation matrix

%% main 
Rmat  = LT1mat*LT1mat';
Rdiag = diag(Rmat);
Cmat = diag(1./sqrt(Rdiag))*Rmat*diag(1./sqrt(Rdiag));

end
