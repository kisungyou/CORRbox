function distance = corr_mmd(corrobj1, corrobj2, varargin)

% CORR_MMD computes maximum mean discrepancy, which is a distance metric on
% the space of probability distributions in the spirit of kernel machines.
% In our implementation, we use the standard gaussian kernel (also known as
% squared exponential kernel) with a limited number of geometries available
% to align with theoretical results for manifold-valued data. 
%
%   * USAGE
%       distance = corr_mmd(input1, input2)
%       distance = corr_mmd(input1, input2, Name, Value)
%
%   * INPUT
%       corrobj1    - an object from 'corr_init' for (p,p,M) data.
%       corrobj2    - an object from 'corr_init' for (p,p,N) data. 
%
%   * PARAMETERS
%       'geometry'  - case-insentitive name of the geometry (default:
%                     "ECM"). Note that it supports a limited number of
%                     options for theoretical reasons. For example, setting
%                     'geometry' with 'qam' will return an error.
%       'bandwidth' - a bandwidth parameter for gaussian kernel. It must be
%                     a positive scalar (default: 1).
%       'biased'    - a logical to determine between biased or unbiased
%                     estimator of the MMD. The unbiased MMD may have
%                     negative values (default: true).
%
%   * OUTPUT
%       distance    - an estimated MMD value.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [10/2023] initial implementation.

%% PREPROCESSING
%  set default values
base_geom  = "ecm";
many_geom  = corrbox_geoms_kernel();
base_scale = 1.0;
base_bias  = true;

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj1", @corraux_checker); 
addRequired(myparser, "corrobj2", @corraux_checker); 

%  parse : parameters
valid_geom  = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
valid_scale = @(x) (isscalar(x) && isnumeric(x) && (x > 0));
addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser,"bandwidth", base_scale, valid_scale);
addParameter(myparser,"biased",base_bias,@islogical);

%  parse
parse(myparser, corrobj1, corrobj2, varargin{:});
parsed = myparser.Results;

%  separate and prepare
input1  = parsed.corrobj1;
input2  = parsed.corrobj2; 

par_geom   = parsed.geometry;
par_sigma  = parsed.bandwidth; 
par_biased = parsed.biased;
geomobj   = corrbox_geomhandle(par_geom);

if (size(input1.data, 1)~=size(input2.data, 1))
    error("* corr_wasserstein : 'corrobj1' and 'corrobj2' should have same dimensions.");
end


%% COMPUTATION
%  pairwise distances
pdmat1 = geomobj.pdist(input1.data);
pdmat2 = geomobj.pdist(input2.data);
dcross = common_pdist2(input1.data, input2.data, par_geom);

%  compute the squared value
dsq = corraux_mmd(pdmat1, pdmat2, dcross, par_sigma, par_biased);

% scale
distance = sqrt(dsq);

end
