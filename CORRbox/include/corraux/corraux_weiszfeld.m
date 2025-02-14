function output = corraux_weiszfeld(C3d, maxiter, thr, Cinit)

%  standard Weiszfeld algorithm on Euclidean space for 3d tensor.
%% computation
%  if one object given, return it
if length(size(C3d)) < 3
    output = C3d;
else
    %  parameters
    P = size(C3d, 1);
    N = size(C3d, 3);

    %  initialization
    if (nargin < 4)
        xold = mean(C3d, 3);
    else
        xold = Cinit;
    end

    %  prep
    dists   = zeros(1,N);
    epsnum  = sqrt(eps);
    weights = ones(1,N)/N; % later can be changed

    %  iteration
    for it=1:maxiter
        % 1. compute the distance
        for n=1:N
            norm2 = norm(xold-C3d(:,:,n),"fro");
            if (norm2 < epsnum)
                dists(n) = norm2 + epsnum;
            else
                dists(n) = norm2;
            end
        end

        % 2. numerator and denominator
        xtmp  = zeros(P,P);
        xtmp2 = 0.0;
        for n=1:N
            xtmp  = xtmp + weights(n)*C3d(:,:,n)/dists(n);
            xtmp2 = xtmp2 + weights(n)/dists(n);
        end
        xnew = xtmp/xtmp2;

        % 3. updating information
        xinc = norm(xnew-xold, "fro");
        xold = xnew;

        % 4. stop/go
        if (xinc < thr)
            break;
        end
    end

    output = xold;
end
end

%   // prepare
%   arma::rowvec xold = xinit;
%   arma::rowvec xtmp(P,fill::zeros);
%   arma::rowvec xnew(P,fill::zeros);
%   arma::vec dists(N,fill::zeros);
%   double xtmp2 = 0.0;
%   double xinc  = 0.0;
%   double norm2 = 0.0;
% 
%   // iteration
%   for (int it=0;it<maxiter;it++){
%     // step 1. compute distance
%     for (int n=0;n<N;n++){
%       norm2 = arma::norm(X.row(n)-xold, 2);
%       if (norm2 < epsnum){
%         dists(n) = norm2 + epsnum;
%       } else {
% 
%     // step 3. updating information
%     xinc = arma::norm(xold-xnew,2);
%     xold = xnew;
%     if (xinc < abstol){
%       break;
%     }
%   }
% 
%   // return
%   return(xold);
% }