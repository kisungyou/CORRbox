function distmat = corr_pdist(corrobj, varargin)

% CORR_PDIST returns a symmetric matrix where each (i,j)-th element is a
% distance measure between i- and j-th elements from given CORR data.
%
%   * USAGE
%       distmat = corr_pdist(corrobj)
%       distmat = corr_pdist(corrobj, Name, Value)
%
%   * INPUT
%       corrobj    - an object from 'corr_init' of size (p,p,N).
%
%   * PARAMETERS 
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%
%   * OUTPUT
%       distmat    - an (N,N) matrix recording distance data pairs.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.
%       0.2. [10/2023] support for multiple geometries added.

%% PREPROCESSING
%  set default values
base_geom = "ecm";
many_geom = corrbox_geoms_all();

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj", @corraux_checker); % using incumbent functions

%  parse : parameters
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
addParameter(myparser, "geometry", base_geom, valid_geom);

%  parse
parse(myparser, corrobj, varargin{:});
parsed = myparser.Results;

%% Computation
%  all inputs
obj3d  = parsed.corrobj;
mygeom = parsed.geometry;

%  derived params
N = size(obj3d.data, 3);
geomobj = corrbox_geomhandle(mygeom);

%  iterate
distmat = zeros(N,N);
for i=1:(N-1)
    tgt1 = obj3d.data(:,:,i);
    for j=(i+1):N
        tgt2 = obj3d.data(:,:,j);
        dval = geomobj.dist(tgt1, tgt2);
        
        distmat(i,j) = dval;
        distmat(j,i) = dval;
    end
end
end