function [label, centers] = corr_kmeans(corrobj, varargin)

% CORR_KMEANS is an adaptation of classical k-means algorithm to
% manifold-valued data on SPD. Lloyd's algorithm was applied in this
% function with an initial configuration obtained from a simple heuristic
% approach for faster computation. 
%
%   * USAGE
%       [label, centers] = corr_kmeans(corrobj)
%       [label, centers] = corr_kmeans(corrobj, Name, Value)
%
%   * INPUT
%       corrobj    - an object from 'corr_init' of size (p,p,N).
%
%   * PARAMETERS 
%       'K'        - predefined number of clusters (default: 2).
%       'maxiter'  - maximum number of iterations to stop (default: 100).
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%
%   * OUTPUT
%       label      - a length-N vector of class label.
%       centers    - (p,p,K) 3d array of cluster centers.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [08/2021] initial implementation.
%       0.2. [10/2023] support for multiple geometries added.

%% PREPROCESSING
%  set default values
base_geom = "ecm";
base_K    = 2;
base_iter = 100;

%  parse : input
myparser = inputParser;
addRequired(myparser, "corrobj", @corraux_checker); % using incumbent functions

%  parse : parameters
many_geom  = corrbox_geoms_all();
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
valid_K    = @(x) (isscalar(x) && isnumeric(x) && (1<x));
valid_iter = @(x) (isscalar(x) && isnumeric(x) && (5<=x));

addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser, "K", base_K, valid_K);
addParameter(myparser, "maxiter", base_iter, valid_iter);

%  parse, separate, and prepare
parse(myparser, corrobj, varargin{:});
parsed = myparser.Results;

myinput = parsed.corrobj; N = size(myinput.data, 3);
myK     = round(parsed.K);
myiter  = max(9, round(parsed.maxiter));
mygeom  = parsed.geometry;

if (myK>=N)
    error("* corr_kmeans : 'K' should be a number in (1,#{data}).");
end

%% COMPUTATION
%  setup the stopping criterion 
par_miter = myiter;
par_incre = 0.01;  

%  set 'old' things
old_label = corraux_kmeans_initialize(myinput, myK);
[old_means, old_variation] = corraux_kmeans_perclass(myinput, old_label, mygeom);
old_cost  = corraux_kmeans_cost(old_variation, old_label);

%  main iteration
itercount = 1;     % record the current iteration count
increment = 10000; % record the increment

while (increment > par_incre)
    %   1. assignment step give old_means
    pdist2    = corraux_pdist2(old_means, myinput.data, mygeom);
    new_label = zeros(size(old_label));
    for i=1:N
        tgtvec = pdist2(:,i);
        tgtidx = find(tgtvec==min(tgtvec));
        if length(tgtidx)>1 % if multiple, choose random one
            new_label(i) = tgtidx(randi(length(tgtidx),1));
        else
            new_label(i) = tgtidx;
        end
    end
    
    %   2. update
    [new_means, new_variation] = corraux_kmeans_perclass(myinput, new_label, mygeom);
    new_cost = corraux_kmeans_cost(new_variation, new_label);
    % Updating is wrong to increase cost function
    % Since it's possible for some runs, only stop it after 25 iterations.
    if ((itercount>=par_miter)&&(new_cost > old_cost)) 
        break;
    end
    increment     = old_cost - new_cost;
    old_means     = new_means;
    old_variation = new_variation;
    old_label     = new_label;
    old_cost      = new_cost;
    
    %   3. stop with itercount
    itercount = itercount + 1;
    if (itercount > par_miter)
        break;
    end
end

%% RESULTS
label  = old_label;
centers = old_means;
% cost = old_cost;

end
