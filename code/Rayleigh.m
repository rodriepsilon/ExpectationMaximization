function y = Rayleigh(x, sigma ) %evaluar la pdf de Rayleigh
    if x<=0
        y = 0; %este valor se obtiene de continuidad en 0
    else
        y = (x/(sigma^2))*exp (-(x^2/(2* sigma^2)));
    end
end