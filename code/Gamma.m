function y = Gamma(x,k,beta) %evaluar la pdf gamma
    if x<=0
        y = 0; %este valor se obtiene de continuidad en 0
    else
      y = ((beta^(k))*(x^(k-1))/gamma(k))*exp(-(beta*x));
    end
end