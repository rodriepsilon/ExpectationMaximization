function y = LogLikelihoodRayleigh(X,Alpha ,Sigma ) 
%log -verosimilitud Rayleigh
    M = length (X);
    N = length (Alpha);
    sum1 = 0;
    for m = 1:M
        sum2 = 0;
        for n = 1:N
            sum2 = sum2 + Alpha (n) * Rayleigh(X(m), Sigma (n));
        end
        sum1 = sum1 + log(sum2);
    end
    y = sum1 ;
end