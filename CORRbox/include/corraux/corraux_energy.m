function Tnm = corraux_energy(DX, DY, DXY)

% https://pages.stat.wisc.edu/~wahba/stat860public/pdf4/Energy/EnergyDistance10.1002-wics.1375.pdf
% parameters
[n,m] = size(DXY);

% we need to compute the following quantities
% compute A : (1/nm)*Sum of all DXY terms
A = (1/(n*m))*sum(sum(DXY));
B = (1/(n*n))*sum(sum(DX));
C = (1/(m*m))*sum(sum(DY));

% compute the statistic
Enm = 2*A - B - C;
Tnm = ((n*m)/(n+m))*Enm;

end