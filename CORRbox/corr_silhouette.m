function score = corr_silhouette(corrobj, label, varargin)

% CORR_SILHOUETTE provides a measure of validating cluster quality defined
% by 'Silhouette' index by Rousseeux (1987). The higher the score is, the 
% better the given clustering is.
%
%   * USAGE
%       score = corr_silhouette(input, label)
%       score = corr_silhouette(input, label, Name, Value)
% 
%   * INPUT
%       corrobj    - an object from 'corr_init' of size (p,p,N).
%       label      - a length-N vector of class labels.
%
%   * PARAMETERS 
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%
%   * OUTPUT
%       score      - a Silhouette score.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation.
%       0.2. [10/2023] support for multiple geometries added.

%% PREPROCESSING
%  set default values
base_geom = "ecm";

%  parse : input
valid_lab = @(x) (isvector(x) && isnumeric(x)); 
myparser  = inputParser;
addRequired(myparser, "corrobj", @corraux_checker); % using incumbent functions
addRequired(myparser, "label", valid_lab);

%  parse : parameters
lists_geom = corrbox_geoms_all();
valid_geom = @(x) (isstring(x) && any(contains(lists_geom, lower(x))));

addParameter(myparser, "geometry", base_geom, valid_geom);

%  parse, separate, and prepare
parse(myparser, corrobj, label, varargin{:});
parsed = myparser.Results;

par_input = parsed.corrobj; N = size(par_input.data, 3);
par_label = round(parsed.label);
par_geom  = parsed.geometry;

if (length(par_label)~=N)
    error('* corr_silhouette : "label" should be a vector of length "corrobj.size(3)".');
end

%% COMPUTATION
%  label control
ulabel = sort(unique(par_label));
K      = length(ulabel);
if (K<2)
    error('* corr_silhouette : the index does not work for K=1.');
end

n       = par_input.size(3);
vec_a   = zeros(n,1); % a(i)
vec_b   = zeros(n,1); % b(i)

indexer = corraux_indexer(par_label);     % index per cluster
distmat = par_geom.pdist(par_input.data); % pairwise distance matrix per geometry
for i=1:n
    % a(i) : average within cluster
    %        singleton index should be taken care.
    idx1     = indexer{(ulabel==par_label(i))};
    if (length(idx1)==1)
        vec_a(i) = 0;
    else
        vec_a(i) = mean(distmat(i,setdiff(idx1,i)));
    end
    % b(i) : weird.. but gotta do.
    otherlabels = setdiff(ulabel, par_label(i)); % non-overlapping classes
    mindistance = zeros(K-1,1);
    for k=1:(K-1)
        idx2 = indexer{(ulabel==otherlabels(k))};
        mindistance(k) = mean(distmat(i,idx2));
    end
    vec_b(i) = min(mindistance);
end

% final computation
vec_s = zeros(n,1);
for i=1:n
    term_top = vec_b(i)-vec_a(i);
    term_bot = max(vec_a(i),vec_b(i));
    vec_s(i) = term_top/term_bot;
end
score = mean(vec_s);

end
