function distmat = corraux_pdist2(array1, array2, geometry)

if nargin < 3
    geometry = "QAM";
end
geomobj = corrbox_geomhandle(geometry);

M = size(array1, 3);
N = size(array2, 3);

distmat = zeros(M,N);
for i=1:M
    tgti = array1(:,:,i);
    for j=1:N
        tgtj = array2(:,:,j);
        distmat(i,j) = geomobj.dist(tgti, tgtj);
    end
end
end