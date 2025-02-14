function coords = ecm_coord_near_mean(C3d)

% given 3d array, find the local coordinates around the mean 

%% computation
%  size
[p,~,N] = size(C3d);

%  compute the frechet mean
[frechet_mean, ~] = ecm_mean(C3d);

%  local coordinates
coords  = zeros(N,round(p*(p-1)/2));
trfmean = trf_COR_to_LT1(frechet_mean);
for n=1:N
    mat_tgt = trf_COR_to_LT1(C3d(:,:,n)) - trfmean;
    coords(n,:) = trf_tril_to_vec(mat_tgt, false);
end

end