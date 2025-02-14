function [Yhat, pred] = corr_kernreg(X, Y, varargin)

% CORR_KERNREG is to fit a kernel regression model for correlation-valued
% data on scalar-valued responses. Given the number of folds, it performs
% cross-validation and hyperparameter optimization for automated learning.
%
%   * USAGE
%       [Yhat, pred] = corr_kernreg(X, Y)
%       [Yhat, pred] = corr_kernreg(X, Y, Name, Value)
%
%   * INPUT
%       X          - an object from 'corr_init' of size (p,p,N).
%       Y          - a length-N vector of response variable.
%
%   * PARAMETERS
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
%       'Kfold'    - the number of folds for cross validation. It should be
%                    a positive integer larger than 1 (default: 5).
%       'parallel' - a logical flag to use parallel computation for
%                    hyperparameter tuning (default: false).
%
%   * OUTPUT
%       Yhat       - a length-N vector of estimated reponses.
%       pred       - a function handle that takes a new object from
%                    'corr_init' consisting of M correlation matrices. If
%                    so, a length-M response vector will be generated as an
%                    output of this function.
%
%   * AUTHOR    Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [10/2023] initial implementation.

%% PREPROCESSING
%  set default values
base_geom = "ecm";
base_fold = 5;
many_geom = corrbox_geoms_all();
base_parallel = false;

%  parse : input
myparser = inputParser;
valid_Y  = @(x) (isvector(x) && isnumeric(x));
addRequired(myparser, "X", @corraux_checker); % using incumbent functions
addRequired(myparser, "Y", valid_Y);         % response should be a vector

%  parse : parameters
valid_geom = @(x) (isstring(x) && any(contains(many_geom, lower(x))));
valid_fold = @(x) (isscalar(x) && isnumeric(x) && (x>1));
valid_parallel = @(x) (isscalar(x) && islogical(x));

addParameter(myparser, "geometry", base_geom, valid_geom);
addParameter(myparser, "Kfold", base_fold, valid_fold);
addParameter(myparser, "parallel", base_parallel, valid_parallel);

%  parse
parse(myparser, X, Y, varargin{:});
parsed = myparser.Results;

par_X = parsed.X; 
par_Y = parsed.Y(:);
if (length(par_Y)~=size(par_X.data,3))
    error("* corr_kernreg : length of 'Y' must be equal to the data size of 'X'.");
end
par_geom = lower(parsed.geometry);
par_fold = max(2, round(parsed.Kfold));
geomobj  = corrbox_geomhandle(par_geom);
par_parallel = parsed.parallel;


%% MAIN
%  compute the Frechet mean
par_iter = 100;
par_thr  = 1e-8;
Cref     = geomobj.fmean(par_X.data, par_iter, par_thr);

%  local coordinates
X_local = geomobj.localise(par_X.data, Cref);

%  train
learned_model = fitrkernel(X_local, par_Y, ...
    'Learner','svm',"OptimizeHyperparameters","auto",...
    'HyperparameterOptimizationOptions',...
    struct('ShowPlots',false,'Verbose',0,"KFold",par_fold,"UseParallel",par_parallel));

%  fit
Yhat = predict(learned_model, X_local);

%  function handle for a new object
pred = @(C3dnew) geomobj.pred_kernreg(C3dnew, Cref, learned_model, "corr_kernreg");

end