function y = ExpectedValueQLognormal (X,Alpha ,Mu ,Sigma , AlphaNew,MuNew , SigmaNew ) %valor esperado de la log - verosimilitud
    N = length ( Alpha );
    M = length (X);
    sum = 0;
    for m = 1:M
        for n = 1:N %el valor esperado Q depende del valor esperado de los datos ocultos
            sum = sum + ExpectedValueYLognormal (n,m,X,Alpha ,Mu, Sigma ) * log ( AlphaNew (n)* lognormal (X(m),MuNew (n),SigmaNew (n)));
        end
    end
    y = sum;
end