function y = ExpectedValueYRayleigh (n,m,X,Alpha, Sigma )
%valor esperado de los datos ocultos
    N = length ( Alpha );
    num = Alpha (n)* Rayleigh(X(m),Sigma (n));
    den = 0;
    for k = 1:N
        den = den + Alpha (k)* Rayleigh(X(m),Sigma(k));
    end
    y = num/den;
end