function output = corr_init(array3d, thr)

% CORR_INIT converts a 3d array into a specific struct that is
% repeatedly used for all exposed functions in our package. If any slice of
% provided 3d array violates one of the requirements - symmetry, unit diagonal, 
% and positive definite, this function tries to adjust the input.
%
%   * USAGE
%       output = corr_init(array3d)
%       output = corr_init(array3d, thr)
%
%   * INPUT
%       array3d     a 3d array of size (p,p,N) of N correlation matrices.
%       thr         threshold value to zero out non-positive eigenvalues,
%
%   * OUTPUT
%       output      a struct containing following elements,
%         .name     name of manifold 'corrbox'
%         .size     a vector of data size
%         .data     (p,p,N) 3d array of processed/checked data
%
%   * AUTHOR     Kisung You (kisungyou@outlook.com)
%   * HISTORY
%       0.1. [07/2021] initial implementation
%       0.2. [10/2023] replace names of the function and returned struct

%% checkers
%   1. 3d array
size3d = size(array3d);
if (length(size3d)~=3)
    error('* corr_init : an input needs to be a 3-dimensional array.');
end
%   2. no need for single slice
if (size3d(3)==1)
    error('* corr_init : single sliced input is not considered.');
end
%   3. square matrices
if (size3d(1)~=size3d(2))
    error('* corr_init : each slice must be squared matrices');
end
p = size3d(1);
N = size3d(3);

%% iterate over slices
if (nargin < 2)
    thr = (p^2)*eps();
end
if (thr < eps())
    error('* corr_init : THR must be a small, positive number.');
end
mydata = zeros(p,p,N);
for i=1:N
    tgt = array3d(:,:,i);
    if (~issymmetric(tgt))
        tgt = (tgt + tgt')/2;
    end
    tmp   = corraux_adjust(tgt, thr);
    mydata(:,:,i) = corraux_cov2cor(tmp);
end

%% return
output = struct();
output.name = 'corrbox';
output.size = [p,p,N];
output.data = mydata;


end