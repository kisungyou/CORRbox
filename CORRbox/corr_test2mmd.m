function pvalue = corr_test2mmd(corrobj1, corrobj2, varargin)

% CORR_TEST2MMD performs a hypothesis testing for the equality of two
% distributions given samples 'corrobj1' and corrobj2' according to the 
% nonparametric method using Maximum Mean Discrepancy.
%
%   * USAGE
%       pvalue = corr_test2mmd(corrobj1, corrobj2)
%       pvalue = corr_test2mmd(corrobj1, corrobj2, Name, Value)
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
%       'bandwidth' - a bandwidth parameter for gaussian kernel. It must be
%                     a positive scalar (default: 1).
%       'biased'    - a logical to determine between biased or unbiased
%                     estimator of the MMD. The unbiased MMD may have
%                     negative values (default: true).
%
%   * OUTPUT
%       pvalue      - p-value under 'H0 : two are equally distributed'.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [11/2023] initial implementation.
%
%   * See also CORR_MMD

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
valid_scale = @(x) (isscalar(x) && isnumeric(x) && (x > 0));

addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser, "niter", base_iter, valid_iter);
addParameter(myparser,"bandwidth", base_scale, valid_scale);
addParameter(myparser,"biased",base_bias,@islogical);

%  parse
parse(myparser, corrobj1, corrobj2, varargin{:});
parsed = myparser.Results;

%  separate and prepare
input1  = parsed.corrobj1; 
input2  = parsed.corrobj2; 
par_iter   = max(9, round(parsed.niter));
par_geom   = parsed.geometry;
par_sigma  = parsed.bandwidth; 
par_biased = parsed.biased;
geomobj    = corrbox_geomhandle(par_geom);

if (size(input1.data, 1)~=size(input2.data, 1))
    error("* corr_test2mmd : 'corrobj1' and 'corrobj2' should have same dimensions.");
end


%% COMPUTATION
%  base distances
DX0 = geomobj.pdist(input1.data); 
DY0 = geomobj.pdist(input2.data);
DZ0 = corraux_pdist2(input1.data, input2.data, mygeom);
DXY = [DX0, DZ0;DZ0', DY0]; % concatenate for future use.

%  test statistic
Tmn = corraux_mmd(DX0, DY0, DZ0, par_sigma, par_biased);

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
    Tvec(i) = corraux_mmd(DX1, DY1, DZ1, par_sigma, par_biased);
end

%% P-VALUE
pvalue = (sum(Tvec>=Tmn)+1)/(par_iter+1);

end