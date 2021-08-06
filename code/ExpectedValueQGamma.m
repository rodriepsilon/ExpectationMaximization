function y = ExpectedValueQGamma(X,Alpha,k,beta,AlphaNew,kNew , betaNew )
%valor esperado de la log - verosimilitud
    N = length ( Alpha );
    M = length (X);
    sum = 0;
    for m = 1:M
        for n = 1:N 
            %el valor esperado Q depende del valor esperado de los datos ocultos
            sum = sum + ExpectedValueYGamma(n,m,X,Alpha ,k,beta ) * log ( AlphaNew (n)* Gamma(X(m),kNew (n),betaNew (n)));
        end
    end
    y = sum;
end