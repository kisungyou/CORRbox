function coords = qam_coord_near_mean(C3d)

% given 3d array, find the local coordinates around the mean 

%% computation
%  size
[p,~,N] = size(C3d);

%  compute the frechet mean
[frechet_mean, ~] = qam_mean(C3d);

%  local coordinates
coords = zeros(N,round(p*(p-1)/2));
for n=1:N
    mat_tgt = corraux_log(frechet_mean, C3d(:,:,n));
    coords(n,:) = trf_tril_to_vec(mat_tgt, false);
end

end