function output = trf_tril_to_vec(L, include_diag)

% given a lower-triangular (LT) matrix, extract LT part and form it as a
% vector for dimension-reduced representation.


% masking of lower-triangular elements
if (include_diag == true)
    % case 1 : include the diagonal part - ECM vectorization
    mask = tril(true(size(L)));
else
    % case 2 : remove diagonal part      - LEC vectorization
    mask = tril(true(size(L)),-1);
end

% Use mask to select lower triangular elements from input array
output = L(mask);

end
