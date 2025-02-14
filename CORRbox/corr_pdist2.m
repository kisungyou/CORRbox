function distmat = corr_pdist2(corrobj1, corrobj2, varargin)

% Given two wrapped data of M and N correlation matrices, CORR_PDIST2
% returns an (M x N) matrix where each (i,j)-th element is a distance
% measure between the i-th observation from the first data and the j-th
% observation from the second data.
%
%   * USAGE
%       distmat = corr_pdist2(corrobj1, corrobj2)
%       distmat = corr_pdist2(corrobj1, corrobj2, Name, Value)
%
%   * INPUT
%       corrobj1   - an object from 'corr_init' for (p,p,M) data.
%       corrobj2   - an object from 'corr_init' for (p,p,N) data. 
%
%   * PARAMETERS
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%
%   * OUTPUT
%       distmat    - an (M,N) matrix recording distance data pairs.
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
addRequired(myparser, "corrobj1", @corraux_checker); 
addRequired(myparser, "corrobj2", @corraux_checker); 

%  parse : parameters
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
addParameter(myparser, "geometry", base_geom, valid_geom);

%  parse
parse(myparser, corrobj1, corrobj2, varargin{:});
parsed = myparser.Results;

%% Computation
%  all inputs
input1   = parsed.corrobj1; 
input2   = parsed.corrobj2;
par_geom = parsed.geometry; 

if (size(input1.data, 1)~=size(input2.data, 1))
    error("* corr_pdist2 : 'corrobj1' and 'corrobj2' should have same dimensions.");
end

%  call
distmat = common_pdist2(input1.data, input2.data, par_geom);

end
