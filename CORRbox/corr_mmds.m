function [embedding, stress] = corr_mmds(corrobj, varargin)

% CORR_MMDS performs a classical multidimensional scaling for visualizing 
% the distribution of N observation on the correlation manifold.
%
%   * USAGE
%       [embedding,stress] = corr_mmds(corrobj)
%       [embedding,stress] = corr_mmds(corrobj, Name, Value)
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
%       stress     - the normalized stress value at the end of iterations.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.
%       0.2. [10/2023] support for multiple geometries added.


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
valid_ndim = @(x) (isscalar(x) && isnumeric(x) && (1<=x) && (x<=min(p,N)));
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));

addParameter(myparser, "ndim", base_dim, valid_ndim);
addParameter(myparser, "geometry", base_geom, valid_geom);

%  parse
parse(myparser, corrobj, varargin{:});
parsed = myparser.Results;


%% COMPUTATION
%  all inputs
data3d = parsed.corrobj;
mygeom = parsed.geometry;
myndim = parsed.ndim;

%  settings
N       = size(data3d.data,3);
geomobj = corrbox_geomhandle(mygeom);

%  compute the pairwise distance matrix
mat_pdist = zeros(N,N);
for i=1:(N-1)
    tgt1 =  data3d.data(:,:,i);
    for j=(i+1):N
        tgt2 = data3d.data(:,:,j);
        dval = geomobj.dist(tgt1, tgt2);
        
        mat_pdist(i,j) = dval;
        mat_pdist(j,i) = dval;
    end
end

% apply 'mdscale' function from Statistics & Machine Learning Toolbox
[embedding, stress] = mdscale(mat_pdist, myndim);


end