function embedding = corr_tsne(corrobj, varargin)

% CORR_TSNE is an extension of the standard t-Stochastic Neighbor Embedding
% algorithm to the Correlation manifold. For numerical efficiency, we use
% the Barnes-Hut algorithm as the default method.
%
%   * USAGE
%       embedding = corr_tsne(corrobj)
%       embedding = corr_tsne(corrobj, Name, Value)
%   
%   * INPUT
%       corrobj       - an object from 'corr_init' of size (p,p,N).
%
%   * PARAMETERS
%       'ndim'        - the target dimension (default: 2).
%       'geometry'    - case-insentitive name of the geometry
%                       (default: "ECM")
%       'perplexity'  - a control parameter to control effective number of
%                       local neighbors (default: 30).
%       'standardize' - a logical to normalize the data (default: false).
%       'exaggerate'  - size of natural clusters in the data (default: 4).
%   
%   * OUTPUT
%       embedding     - an (N,ndim) embedding coordinates.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.
%       0.2. [10/2023] support for multiple geometries added.
%                      simplified its output to only deliver embeddings.

%% PREPROCESSING
%  set default values
base_dim  = 2;
base_geom = "ecm";
many_geom = corrbox_geoms_all();

%  special parameters for tsne that would be passed
base_perplexity = 30;     % positive scalar
base_standardize = false; % boolean
base_exaggerate  = 4;     % scalar value of 1 or greater

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj", @corraux_checker); % using incumbent functions
[p,~,N] = size(corrobj.data);

%  parse : parameters
valid_ndim = @(x) (isscalar(x) && isnumeric(x) && (1<=x) && (x<=min(p,N)));
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
valid_perplexity = @(x) (isscalar(x) && isnumeric(x) && (x>0));
valid_exaggerate = @(x) (isscalar(x) && isnumeric(x) && (x>=1));

addParameter(myparser, "ndim", base_dim, valid_ndim);
addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser, "perplexity", base_perplexity, valid_perplexity);
addParameter(myparser, "standardize", base_standardize, @islogical);
addParameter(myparser, "exaggerate", base_exaggerate, valid_exaggerate);

%  parse, separate, and prepare
parse(myparser, corrobj, varargin{:});
parsed = myparser.Results;

myinput = parsed.corrobj; 
mygeom  = parsed.geometry;
geomobj = corrbox_geomhandle(mygeom);

par_ndim  = max(1, round(parsed.ndim));
par_perp  = parsed.perplexity;
par_scale = parsed.standardize;
par_exagg = parsed.exaggerate;


%% MAIN
%  compute the local coordinates
X = geomobj.localcoord(myinput.data);

%  run through the matlab routine
embedding = tsne(X, 'NumDimensions', par_ndim, ...
    'Perplexity',par_perp,'Standardize',par_scale,...
    'Exaggeration',par_exagg);

end