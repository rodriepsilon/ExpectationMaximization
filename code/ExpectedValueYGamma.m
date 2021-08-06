function y = ExpectedValueYGamma(n,m,X,Alpha,k,beta) 
%valor esperado de los datos invisibles
    N = length ( Alpha );
    num = Alpha (n)*Gamma(X(m),k(n),beta(n));
    den = 0;
    for i = 1:N
        den = den + Alpha(i)* Gamma(X(m),k(i),beta(i));
    end
    y = num/den;
end