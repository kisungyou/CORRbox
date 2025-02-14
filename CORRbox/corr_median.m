function [center, variation] = corr_median(corrobj, varargin)

% CORR_MEDIAN computes Fréchet median and its variation given an object
% processed by 'corr_init' function. Variation is analogous to the sample
% version of mean absolute deviation for the data on CORR manifold.
%
%   * USAGE
%       [center, variation] = corr_median(input)
%       [center, variation] = corr_median(input, Name, Value)
%
%   * INPUT
%       corrobj    - an object from 'corr_init' of size (p,p,N).
%
%   * PARAMETERS
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%       'maxiter'  - the number of maximum iterations (default: 100).
%       'thr'      - stopping criterion for iterations (default: 1e-8).
%
%   * OUTPUT
%       center     - an empirical Fréchet median matrix of size (p,p).
%       variation  - sample Fréchet variation.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.
%       0.2. [10/2023] support for multiple geometries added.

%% PREPROCESSING
%  set default values
base_geom = "ecm";
many_geom = corrbox_geoms_all();
base_iter = 100;
base_thr  = 1e-8;

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj", @corraux_checker); % using incumbent functions

%  parse : parameters
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
valid_iter = @(x) (isscalar(x) && isnumeric(x) && (2 < x));
valid_thr  = @(x) (isscalar(x) && isnumeric(x) && (100*eps < x));

addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser, "maxiter",  base_iter, valid_iter);
addParameter(myparser, "thr", base_thr, valid_thr);

%  parse
parse(myparser, corrobj, varargin{:});
parsed = myparser.Results;

%% COMPUTATION
%  all inputs
obj3d  = parsed.corrobj;
mygeom = parsed.geometry;
myiter = parsed.maxiter;
mythr  = parsed.thr;

%  settings
geomobj = corrbox_geomhandle(mygeom);

%  main calling
[center, variation] = geomobj.fmedian(obj3d.data, myiter, mythr);

end