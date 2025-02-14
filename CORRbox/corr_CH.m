function score = corr_CH(corrobj, label, varargin)


% corr_CH is a measure of clustering quality that was proposed by 
% Calinski and Harabasz in 1974. We interpret the score in a way that the
% higer the score is, the better the given clustering is.
%
%   * USAGE
%       score = corr_CH(input, label)
%       score = corr_CH(input, label, Name, Value)
% 
%   * INPUT
%       corrobj    - an object from 'corr_init' of size (p,p,N).
%       label      - a length-N vector of class labels.
%
%   * PARAMETERS 
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%
%   * OUTPUT
%       score      - a CH score.
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

par_input = parsed.corrobj; 
N = size(par_input.data, 3);
p = size(par_input.data, 1);
par_label = round(parsed.label);
par_geom  = parsed.geometry;

if (length(par_label)~=N)
    error('* corr_CH : "label" should be a vector of length "corrobj.size(3)".');
end


%% MAIN COMPUTATION
%  Preliminary
ulabel = sort(unique(label));
K      = length(ulabel);
if (K<2)
    error('* corr_CH : index does not work for K=1.');
end

%  Compute Means
mean_group = zeros(p,p,K);
for k=1:K
    idk = find(label==ulabel(k));
    [mean_tmp,~] = par_geom.fmean(par_input.data(:,:,idk));
    mean_group(:,:,k) = mean_tmp;
end
[mean_all,~]   = par_geom.fmean(par_input.data);

% Compute Separation (numerator)
val_numerator = 0;
for k=1:K
    val_numerator = val_numerator + (par_geom.dist(mean_all, mean_group(:,:,k))^2)*sum(label==ulabel(k))/(K-1);
end

% Compute Cohesion   (denominator)
val_denominator = 0;
for k=1:K
    tmp_sum  = 0;
    tmp_id   = find(label==ulabel(k));
    tmp_nobs = length(tmp_id);
    for it=1:tmp_nobs
        tmp_sum = tmp_sum + (par_geom.dist(mean_group(:,:,k), input.data(:,:,tmp_id(it)))^2)/(N-K);
    end
    val_denominator = val_denominator + tmp_sum;
end

% return the score 
score = val_numerator/val_denominator;


end
