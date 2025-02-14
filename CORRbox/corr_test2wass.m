function pvalue = corr_test2wass(corrobj1, corrobj2, varargin)

% CORR_TEST2WASS performs a hypothesis testing for the equality of two
% distributions given samples 'corrobj1' and corrobj2' according to the 
% nonparametric method using 2-Wasserstein distance.
%
%   * USAGE
%       pvalue = corr_test2wass(corrobj1, corrobj2)
%       pvalue = corr_test2wass(corrobj1, corrobj2, Name, Value)
%
%   * INPUT
%       corrobj1   - an object from 'corr_init' for (p,p,M) data.
%       corrobj2   - an object from 'corr_init' for (p,p,N) data. 
%
%   * PARAMETERS
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%       'niter'    - the number of iterations (default: 999).
%
%   * OUTPUT
%       pvalue     - p-value under 'H0 : two are equally distributed'.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.
%       0.2. [10/2021] fixed a critical error.
%       0.3. [10/2023] support for multiple geometries added.


%% PREPROCESSING
%  set default values
base_geom = "ecm";
many_geom = corrbox_geoms_all();
base_iter  = 999;

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj1", @corraux_checker); 
addRequired(myparser, "corrobj2", @corraux_checker); 

%  parse : parameters
valid_geom  = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
valid_iter  = @(x) (isnumeric(x) && isscalar(x) && (9<=x));

addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser, "niter", base_iter, valid_iter);

%  parse
parse(myparser, corrobj1, corrobj2, varargin{:});
parsed = myparser.Results;

%  separate and prepare
input1  = parsed.corrobj1; m = size(input1.data, 3);
input2  = parsed.corrobj2; n = size(input2.data, 3);
myiter  = max(9, round(parsed.niter));
mygeom  = parsed.geometry;
geomobj = corrbox_geomhandle(mygeom);

if (size(input1.data, 1)~=size(input2.data, 1))
    error("* corr_test2wass : 'corrobj1' and 'corrobj2' should have same dimensions.");
end


%% COMPUTATION
%  base distances
DX0 = geomobj.pdist(input1.data); 
DY0 = geomobj.pdist(input2.data);
DZ0 = corraux_pdist2(input1.data, input2.data, mygeom);
DXY = [DX0, DZ0;DZ0', DY0]; % concatenate for future use.

%  test statistic
Tmn = corraux_wassersteinD(DZ0, 2.0);

%  main iterations
Tvec = zeros(1,myiter);
for i=1:myiter
    % random permutation
    id0 = randperm(m+n); % random permutation
    idx = id0(1:m);
    idy = id0((m+1):(m+n));
    
    DX1 = DXY(idx,idx);
    DY1 = DXY(idy,idy);
    DZ1 = DXY(idx,idy);
    Tvec(i) = corraux_wassersteinD(DZ1, 2.0);
end

%% P-VALUE
pvalue = (sum(Tvec>=Tmn)+1)/(myiter+1);

end