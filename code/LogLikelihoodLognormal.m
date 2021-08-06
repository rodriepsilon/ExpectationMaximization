function y = LogLikelihoodLognormal (X,Alpha ,Mu , Sigma ) %log -verosimilitud
    M = length (X);
    N = length ( Alpha );
    sum1 = 0;
    for m = 1:M
        sum2 = 0;
        for n = 1:N
            sum2 = sum2 + Alpha (n) * lognormal (X(m),Mu(n),Sigma (n));
        end
        sum1 = sum1 + log ( sum2 );
    end
    y = sum1 ;
end