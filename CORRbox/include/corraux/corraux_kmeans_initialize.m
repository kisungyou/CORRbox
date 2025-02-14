function initlabel = corraux_kmeans_initialize(input, K)

% initialize the label for any iterative clustering algorithm using a
% simple, proxy geometric operations combining pseudo-equivariant embedding
% and random projection.

% size
p = input.size(1);
N = input.size(3);

% data reformulation
mydata = zeros(N, p*p);
for i=1:N
    tgtlg       = logm(input.data(:,:,i)); % equivariant embedding
    mydata(i,:) = reshape(tgtlg, [1, p*p]);
end

% transform into low-dimensional space if original dimensionality is higher
% than some threshold; this will enhance stability. 
if size(mydata, 2) > 5
    mydata = mydata*randn(p*p, 5);
end

% apply matlab's original kmeans function
initlabel = kmeans(mydata,K);

end