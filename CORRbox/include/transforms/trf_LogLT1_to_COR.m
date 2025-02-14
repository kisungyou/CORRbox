function Cmat = trf_LogLT1_to_COR(LT1log)

% input
%   LT1log : lower-triangular with unit diagonals with logm applied
% output
%   Cmat   : correlation matrix

%% main 
LT1mat = expm(LT1log);
Rmat   = LT1mat*LT1mat';
Rdiag  = diag(Rmat);
Cmat   = diag(1./sqrt(Rdiag))*Rmat*diag(1./sqrt(Rdiag));

end
