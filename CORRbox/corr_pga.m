function [embedding,pred] = corr_pga(corrobj, varargin)

% CORR_PGA is an implementation of the method of (tangential) principal
% geodesic analysis, an extension of the celebrated PCA to the
% manifold-valued data. 
%
%   * USAGE
%       [embedding, pred] = corr_pga(corrobj)
%       [embedding, pred] = corr_pga(corrobj, Name, Value)
%   
%   * INPUT
%       corrobj    - an object from 'corr_init' of size (p,p,N).
%
%   * PARAMETERS
%       'ndim'     - the target dimension (default: 2).
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%   
%   * OUTPUT
%       embedding  - an (N,ndim) embedding coordinates.
%       pred       - a function handle that takes a new object from
%                    'corr_init' consisting of M correlation matrices. As
%                    an output of the pred function, an (M,ndim) matrix of
%                    embeddings for the new input is returned.
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

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj", @corraux_checker); % using incumbent functions
[p,~,N] = size(corrobj.data);

%  parse : parameters
valid_ndim = @(x) (isscalar(x) && (1<=x) && (x<=min(p,N)));
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));

addParameter(myparser, "ndim", base_dim, valid_ndim);
addParameter(myparser, "geometry", base_geom, valid_geom);

%  parse, separate, and prepare
parse(myparser, corrobj, varargin{:});
parsed = myparser.Results;

myinput = parsed.corrobj; 
myndim  = max(1, round(parsed.ndim));
mygeom  = parsed.geometry;

%% COMPUTATION
%  main embedding routine
[embedding, Cmean, apply_mean, apply_proj] = common_pga_main(myinput.data, myndim, mygeom);

%  prediction
pred = @(newobj) common_pga_pred(newobj, Cmean, apply_mean, apply_proj, mygeom);

end