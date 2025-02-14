function embedding = qam_pga(C3d, ndim)

%  compute the intrinsic mean with arbitrary control parameters 
[p,~,N]  = size(C3d);
par_iter = 100;
par_thr  = 1e-6;
[center, ~] = qam_mean(C3d, par_iter, par_thr);

%  find the local coordinates
stretched = zeros(N,p^2); % tangentialize
for i=1:N
    logged = real(corraux_log(center, C3d(:,:,i)));
    logged = (logged + logged')/2;
    stretched(i,:) = logged(:); % stack as rows
end

%  apply PCA near the mean
embedding = corraux_PCA(stretched, round(ndim));

end