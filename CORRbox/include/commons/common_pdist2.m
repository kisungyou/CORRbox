function crossdist = common_pdist2(X3, Y3, geometry)

% parameters
M = size(X3, 3);
N = size(Y3, 3);
geomobj = corrbox_geomhandle(geometry);

% compute
crossdist = zeros(M,N);
for m=1:M
    tgt1 = X3(:,:,m);
    for n=1:N
        tgt2 = Y3(:,:,n);
        crossdist(m,n) = geomobj.dist(tgt1, tgt2);
    end
end
end