function [embedding,pred] = corr_ae(corrobj, varargin)

% CORR_AE incorporates a simple autoencoder with a single hidden layer
% under the extrinsic regime. The current version incorporates two schemes 
% of L2-weight and sparse regularization for robust learning of latent
% representation for correlation-valued data. 
%
%   * USAGE
%       [embedding, pred] = corr_ae(corrobj)
%       [embedding, pred] = corr_ae(corrobj, Name, Value)
%
%   * INPUT
%       corrobj    - an object from 'corr_init' of size (p,p,N).
%
%   * PARAMETERS
%       'ndim'     - the target dimension (default: 2).
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%       'alpha'    - a parameter that controls L2-weight regularization
%                    (default: 0.01). It should be larger than 0. 
%       'beta'     - a parameter that controls sparse regularization of
%                    network's weights (default: 1). Must be larger than 0.
%
%   * OUTPUT
%       embedding  - an (N,ndim) embedding coordinates.
%       pred       - a function handle that takes a new object from
%                    'corr_init' consisting of M correlation matrices. As
%                    an output of the pred function, an (N,ndim) matrix of
%                    embeddings for the new input is returned.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [10/2023] initial implementation.

%% PREPROCESSING
%  set default values
base_dim  = 2;
base_geom = "ecm";
many_geom = corrbox_geoms_all();
base_alpha = 0.01;
base_beta  = 1;

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj", @corraux_checker); % using incumbent functions
[p,~,N] = size(corrobj.data);

%  parse : parameters
valid_ndim = @(x) (isscalar(x) && (1<=x) && (x<=min(p,N)));
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
valid_lambda = @(x) (isscalar(x) && isnumeric(x) && (x>0));

addParameter(myparser, "ndim", base_dim, valid_ndim);
addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser, "alpha", base_alpha, valid_lambda);
addParameter(myparser, "beta", base_beta, valid_lambda);

%  parse, separate, and prepare
parse(myparser, corrobj, varargin{:});
parsed = myparser.Results;

myinput = parsed.corrobj; 
myndim  = max(1, round(parsed.ndim)); % this will correspond to hidden 
mygeom  = parsed.geometry;
myalpha = parsed.alpha;
mybeta  = parsed.beta;

%% MAIN
%  call the main routine - no need for separating
[embedding, pred] = common_ae_main(myinput.data, myndim, mygeom, myalpha, mybeta);


end