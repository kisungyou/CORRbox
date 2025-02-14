function distance = corr_wasserstein(corrobj1, corrobj2, varargin)

% CORR_WASSERSTEIN aims at computing the Wasserstein distance of order 'p'
% between two empirical measures of correlation matrices. 
%
%   * USAGE
%       distance = corr_wasserstein(input1, input2)
%       distance = corr_wasserstein(input1, input2, Name, Value)
%
%   * INPUT
%       corrobj1   - an object from 'corr_init' for (p,p,M) data.
%       corrobj2   - an object from 'corr_init' for (p,p,N) data. 
%
%   * PARAMETERS
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%       'order'    - the order of Wasserstein metric (default: 2).
%
%   * OUTPUT
%       distance   - an estimated p-Wasserstein distance.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [10/2023] initial implementation.
%       0.2. [10/2023] code refactored to match formal MATLAB styles.

%% PREPROCESSING
%  set default values
base_geom = "ecm";
many_geom = corrbox_geoms_all();
base_order = 2;

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj1", @corraux_checker); 
addRequired(myparser, "corrobj2", @corraux_checker); 

%  parse : parameters
valid_geom  = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
valid_order = @(x) (isscalar(x) && isnumeric(x) && (1<=x));
addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser, "order", base_order, valid_order);

%  parse
parse(myparser, corrobj1, corrobj2, varargin{:});
parsed = myparser.Results;

%  separate and prepare
input1  = parsed.corrobj1; M = size(input1.data, 3);
input2  = parsed.corrobj2; N = size(input2.data, 3);
mygeom  = parsed.geometry;
myp     = double(parsed.order);
geomobj = corrbox_geomhandle(mygeom);

if (size(input1.data, 1)~=size(input2.data, 1))
    error("* corr_wasserstein : 'corrobj1' and 'corrobj2' should have same dimensions.");
end


%% COMPUTATION
%  pairwise distance computation
distmat = zeros(M,N);
for m=1:M
    for n=1:N
        distmat(m,n) = geomobj.dist(input1.data(:,:,m), input2.data(:,:,n));
    end
end

%  create weight vectors (currently, uniform)
wx = ones(M,1)/M; wx = wx(:);
wy = ones(N,1)/N; wy = wy(:);

% invoke a subroutine
distance = corraux_wassersteinD(distmat, myp, wx, wy);

end
