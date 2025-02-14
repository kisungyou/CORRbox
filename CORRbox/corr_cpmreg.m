function [Yhat, pred] = corr_cpmreg(X, Y, varargin)

% CORR_CPMREG is a rusty implementation of connectome-based predictive
% modeling without cross validation.
%
%   * USAGE
%       [Yhat, pred] = corr_cpmreg(X, Y)
%
%   * INPUT
%       X          - an object from 'corr_init' of size (p,p,N).
%       Y          - a length-N vector of response variable.
%                    hyperparameter tuning (default: false).
%   * PARAMETERS
%       'geometry' - case-insentitive name of the geometry (default: "ECM"). 
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
base_thr  = 0.05;
valid_thr = @(x) ((x<=1) && (x>=0));

%  parse : input
myparser = inputParser;
valid_Y  = @(x) (isvector(x) && isnumeric(x));
addRequired(myparser, "X", @corraux_checker); % using incumbent functions
addRequired(myparser, "Y", valid_Y);         % response should be a vector

%  parse : parameters
addParameter(myparser, "threshold", base_thr, valid_thr);


%  parse
parse(myparser, X, Y);
parsed = myparser.Results;

par_X = parsed.X; 
par_Y = parsed.Y(:);
if (length(par_Y)~=size(par_X.data,3))
    error("* corr_cpmreg : length of 'Y' must be equal to the data size of 'X'.");
end
par_thr = parsed.threshold;

%% TRAINING
%  get the size
[P,~,N] = size(par_X.data);

%  compute the mask
idxPval = [];
for i=1:(P-1)
    for j=(i+1):P
        % select the data
        now_data = par_X.data(i,j,:);
        now_data = now_data(:);

        % compute the p-value
        [~, now_pval] = corrcoef(now_data, par_Y);

        % conditionally accept it
        if (now_pval(1,2) < par_thr)
            idxPval = [idxPval; i, j];
        end
    end
end

%  subset the data
n_selected = size(idxPval, 1);

%  rearrange the inputs
X_arranged = zeros(N, n_selected);
for n=1:N
    for m=1:n_selected
        X_arranged(n,m) = par_X.data(idxPval(m,1),idxPval(m,2),n);
    end
end

%  fit linear regression model
mdl_lr = fitlm(X_arranged, par_Y);

%  return a fitted value
Yhat = mdl_lr.Fitted;

%% Testing
%  create a prediction function 

    function new_output = cpmreg_predict(new_input, mdl_train, idx_sel)
    % check the input
    if (~corraux_checker(new_input))
        error("* corr_cpmreg : new input is not a valid object.");
    end
    
    % get the size
    [~, ~, new_N] = size(new_input.data);
    new_nsel = size(idx_sel, 1);

    % rearrange the input
    new_arranged = zeros(new_N, new_nsel);
    for new_n=1:new_N
        for new_m=1:new_nsel
            new_arranged(new_n,new_m) = new_input.data(idx_sel(new_m,1), idx_sel(new_m,2), new_n);
        end
    end

    % predict
    new_output = predict(mdl_train, new_arranged);
    end

% create a handle
pred = @(C3dnew) cpmreg_predict(C3dnew, mdl_lr, idxPval);

end