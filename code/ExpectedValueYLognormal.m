function y = ExpectedValueYLognormal (n,m,X,Alpha ,Mu , Sigma ) %valor esperado de los datos invisibles
    N = length ( Alpha );
    num = Alpha (n)* lognormal (X(m),Mu(n),Sigma (n));
    den = 0;
    for k = 1:N
        den = den + Alpha (k)* lognormal (X(m),Mu(k),Sigma (k));
    end
    y = num/den;
end