function y = ExpectedValueQRayleigh(X,Alpha,Sigma , AlphaNew,SigmaNew ) 
%valor esperado de la log - verosimilitud
    N = length ( Alpha );
    M = length (X);
    sum = 0;
    for m = 1:M
        for n = 1:N %el valor esperado Q depende del valor esperado de los datos ocultos
            sum = sum + ExpectedValueYRayleigh(n,m,X,Alpha, Sigma ) * log ( AlphaNew (n)* Rayleigh(X(m),SigmaNew (n)));
        end
    end
    y = sum;
end