function LT1vec = aux_chol2vec(LT1mat)

% convert a LT1 matrix into a vectorial form
mask   = tril(true(size(LT1mat)));
LT1vec = LT1mat(mask);

end