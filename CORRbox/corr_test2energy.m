function pvalue = corr_test2energy(corrobj1, corrobj2, varargin)

% CORR_TEST2ENERGY performs a hypothesis testing for the equality of two
% distributions given samples 'corrobj1' and corrobj2' according to the 
% nonparametric method using Energy distance.
%
%   * USAGE
%       pvalue = corr_test2energy(corrobj1, corrobj2)
%       pvalue = corr_test2energy(corrobj1, corrobj2, Name, Value)
%
%   * INPUT
%       corrobj1    - an object from 'corr_init' for (p,p,M) data.
%       corrobj2    - an object from 'corr_init' for (p,p,N) data. 
%
%   * PARAMETERS
%       'geometry'  - case-insentitive name of the geometry. Note that it
%                     supports a limited number of options for theoretical
%                     reasons. For example, setting 'geometry' as 'QAM'
%                     will return an error (default: "ECM").
%       'niter'     - the number of iterations (default: 999).
%
%   * OUTPUT
%       pvalue      - p-value under 'H0 : two are equally distributed'.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [03/2024] initial implementation.

%% PREPROCESSING
%  set default values
base_geom = "ecm";
many_geom = corrbox_geoms_kernel();

base_iter  = 999;
base_scale = 1.0;
base_bias  = true;

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
input1  = parsed.corrobj1; 
input2  = parsed.corrobj2; 
par_iter   = max(9, round(parsed.niter));
par_geom   = parsed.geometry;
geomobj    = corrbox_geomhandle(par_geom);

if (size(input1.data, 1)~=size(input2.data, 1))
    error("* corr_test2energy : 'corrobj1' and 'corrobj2' should have same dimensions.");
end


%% COMPUTATION
%  base distances
DX0 = geomobj.pdist(input1.data); 
DY0 = geomobj.pdist(input2.data);
DZ0 = corraux_pdist2(input1.data, input2.data, mygeom);
DXY = [DX0, DZ0;DZ0', DY0]; % concatenate for future use.

%  test statistic
Tmn = corraux_energy(DX0, DY0, DZ0);

%  main iterations
Tvec = zeros(1,par_iter);
for i=1:par_iter
    % random permutation
    id0 = randperm(m+n); % random permutation
    idx = id0(1:m);
    idy = id0((m+1):(m+n));
    
    DX1 = DXY(idx,idx);
    DY1 = DXY(idy,idy);
    DZ1 = DXY(idx,idy);
    Tvec(i) = corraux_energy(DX1, DY1, DZ1);
end

%% P-VALUE
pvalue = (sum(Tvec>=Tmn)+1)/(par_iter+1);

end