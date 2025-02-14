function [label, centers] = corr_kmedoids(corrobj, varargin)

% CORR_KMEDOIDS performs a k-medoids clustering on the given sample of 
% N observations on the correlation manifold.
%
%   * USAGE
%       [label, centers] = corr_kmedoids(corrobj)
%       [label, centers] = corr_kmedoids(corrobj, Name, Value)
%
%   * INPUT
%       corrobj    - an object from 'corr_init' of size (p,p,N).
%
%   * PARAMETERS 
%       'K'        - predefined number of clusters (default: 2).
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%
%   * OUTPUT
%       label      - a length-N vector of class label.
%       centers    - (p,p,K) 3d array of cluster centers.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.
%       0.2. [10/2023] support for multiple geometries added.

%% PREPROCESSING
%  set default values
base_geom = "ecm";
many_geom = corrbox_geoms_all();
base_K    = 2;

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj", @corraux_checker); % using incumbent functions

%  parse : parameters
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
valid_K    = @(x) (isscalar(x) && isnumeric(x) && (1<x));
addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser, "K", base_K, valid_K);

%  parse, separate, and prepare
parse(myparser, corrobj, varargin{:});
parsed = myparser.Results;

myinput = parsed.corrobj; N = size(myinput.data, 3);
myK     = round(parsed.K);
mygeom  = parsed.geometry;
geomobj = corrbox_geomhandle(mygeom);

if (myK>=N)
    error("* corr_kmedoids : 'K' should be a number in (1,#{data}).");
end

%% COMPUTATION
%  Since squareform is not supported, we need local coordinates.
new_data = geomobj.localcoord(myinput.data);

%  apply k-medoids
[label,~,~,~,midx] = kmedoids(new_data, myK);

%  extract centroids from PAM
centers = myinput.data(:,:,midx);


end
