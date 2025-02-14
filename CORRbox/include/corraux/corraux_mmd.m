function dsq = corraux_mmd(pdmat1, pdmat2, crossdmat, par_sigma, par_biased)

% return a squared value for the bias handling

% apply Gaussian kernel
exp1  = exp(-(pdmat1.^2)/(2*(par_sigma^2)));
exp2  = exp(-(pdmat2.^2)/(2*(par_sigma^2)));
exp12 = exp(-(crossdmat.^2)/(2*(par_sigma^2)));

[M,N] = size(exp12);

% case branching
if (par_biased)
    % case 1 : biased
    term1 = sum(exp1,"all")/(M^2);
    term2 = sum(exp2,"all")/(N^2);
    term3 = (2/(M*N))*sum(exp12,"all");
    dsq   = term1+term2-term3;
else
    % case 2 : unbiased
    % prelim to fill zeros in the diagonal
    for m=1:M
        exp1(m,m) = 0;
    end
    for n=1:N
        exp2(n,n) = 0;
    end
    % compute terms
    term1 = sum(exp1,"all")/(M*(M-1));
    term2 = sum(exp2,"all")/(N*(N-1));
    term3 = (2/(M*N))*sum(exp12,"all");
    dsq   = term1+term2-term3;
end
end