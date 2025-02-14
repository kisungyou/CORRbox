function label = corr_specc(corrobj, varargin)

% CORR_SPECC applies a spectral clustering algorithm by a version of 
% Zelnik-Manor and Perona (2005) where an affinity matrix is constructed 
% using the local distance to the 'nnbd'-th nearest observations for N 
% observations on the correlation manifold.
%
%   * USAGE
%       label = corr_specc(input)
%       label = corr_specc(input, Name, Value)
%
% 
%   * INPUT
%       corrobj    - an object from 'corr_init' of size (p,p,N).
%
%   * PARAMETERS 
%       'K'        - predefined number of clusters (default: 2).
%       'nnbd'     - size of the neighborhood for bandwidth (default: 5).
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%
%   * OUTPUT
%       label      - a length-N vector of class label.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [08/2021] initial implementation.
%       0.2. [10/2023] support for multiple geometries added.

%% PREPROCESSING
%  set default values
base_geom = "ecm";
many_geom = corrbox_geoms_all();
base_K    = 2;
base_nnbd = 5;

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj", @corraux_checker); % using incumbent functions

%  parse : parameters
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
valid_K    = @(x) (isnumeric(x) && isscalar(x) && (1<x));
valid_nnbd = @(x) (isnumeric(x) && isscalar(x) && (1<=x));

addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser, "K", base_K, valid_K);
addParameter(myparser, "nnbd", base_nnbd, valid_nnbd);

%  parse, separate, and prepare
parse(myparser, corrobj, varargin{:});
parsed = myparser.Results;

myinput = parsed.corrobj; N = size(myinput.data, 3);
myK     = max(2, round(parsed.K));
mygeom  = parsed.geometry;
mynnbd  = max(1, round(parsed.nnbd));
geomobj = corrbox_geomhandle(mygeom);

if (myK>=N)
    error("* corr_specc : 'K' should be a number in (1,#{data}).");
end
if (mynnbd >= (N-1))
    error("* corr_specc : 'nnbd' is too large. Make it smaller.");
elseif (mynnbd < 2)
    error("* corr_specc : 'nnbd' is too small. Make it larger");
end

%% COMPUTATION
%  compute the pairwise distance matrix
mat_pdist = zeros(N,N);
for i=1:(N-1)
    tgt1 = myinput.data(:,:,i);
    for j=(i+1):N
        tgt2 = myinput.data(:,:,j);
        dval = geomobj.dist(tgt1, tgt2);
        
        mat_pdist(i,j) = dval;
        mat_pdist(j,i) = dval;
    end
end


%  compute the nearest distance vector
vec_nearest = zeros(1,N);
for n=1:N
    tgt_sorted     = sort(mat_pdist(n,:));
    vec_nearest(n) = tgt_sorted(mynnbd);
end

%  compute the similarity matrix using Zelnik-Manor
mat_S = ones(N,N);
for i=1:(N-1)
    for j=(i+1):N
        mat_S(i,j) = exp(-(mat_pdist(i,j)^2)/(vec_nearest(i)*vec_nearest(j)));
        mat_S(j,i) = mat_S(i,j);
    end
end

%  perform spectral clustering using NJW
label = spectralcluster(mat_S, myK, 'Distance','precomputed','LaplacianNormalization','symmetric');


end