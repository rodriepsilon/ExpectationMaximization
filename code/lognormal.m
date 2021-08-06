function y = lognormal (x,mu , sigma ) %evaluar la pdf log - normal
    if x<=0
        y = 0; %este valor se obtiene de continuidad en 0
    else
        y = 1/(x* sigma * sqrt (2* pi))*exp (-( log(x)-mu) ^2/(2* sigma^2));
    end
end