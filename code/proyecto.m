g = table; %tabla que almacena resultados de mezcla de Rayleigh
r = table; %tabla que almacena resultados de mezcla de Gamma
h = table; %tabla que almacena resultados de mezcla de Log-normal

datos = readtable('../data/tidy_data/datos.csv');
for i = 1:36
    display(2*i);
    u = datos{:,2*i};
    v = RayleighMix.fit(u,3);
    column1 = datos.Properties.VariableNames(2*i);
    column2 = v.LogLikelihood;
    column3 = v.ComponentProportions(1);
    column4 = v.ComponentProportions(2);
    column5 = v.ComponentProportions(3);
    column6 = v.sigma(1);
    column7 = v.sigma(2);
    column8 = v.sigma(3);

    tempt = table(column1, column2,column3,column4,column5,column6,column7,column8, 'VariableNames',{'Variable','Logverosimilitud','alfa1', 'alfa2','alfa3','sigma1', 'sigma2','sigma3'});
    r = [r;tempt];
end
writetable(r,'../data/tidy_data/rayleigh.csv')

for i = 1:36
    display(2*i);
    u = datos{:,2*i};
    v = gammaMix.fit(u,3);
    column1 = datos.Properties.VariableNames(2*i);
    column2 = v.LogLikelihood;
    column3 = v.ComponentProportions(1);
    column4 = v.ComponentProportions(2);
    column5 = v.ComponentProportions(3);
    column6 = v.k(1);
    column7 = v.k(2);
    column8 = v.k(3);
    column9 = v.beta(1);
    column10 = v.beta(2);
    column11 = v.beta(3);
    tempt = table(column1, column2,column3,column4,column5,column6,column7,column8,column9,column10,column11, 'VariableNames',{'Variable','Logverosimilitud','alfa1', 'alfa2','alfa3','k1', 'k2','k3','beta1', 'beta2', 'beta3'});
    g = [g;tempt];
end
writetable(g,'../data/tidy_data/gamma.csv')
for i = 1:36
    display(2*i);
    u = datos{:,2*i};
    v = lognormalMix.fit(u,3);
    column1 = datos.Properties.VariableNames(2*i);
    column2 = v.LogLikelihood;
    column3 = v.ComponentProportions(1);
    column4 = v.ComponentProportions(2);
    column5 = v.ComponentProportions(3);
    column6 = v.mu(1);
    column7 = v.mu(2);
    column8 = v.mu(3);
    column9 = v.sigma(1);
    column10 = v.sigma(2);
    column11 = v.sigma(3);
    tempt = table(column1, column2,column3,column4,column5,column6,column7,column8,column9,column10,column11, 'VariableNames',{'Variable','Logverosimilitud','alfa1', 'alfa2','alfa3','mu1', 'mu2','mu3','sigma1', 'sigma2', 'sigma3'});
    h = [h;tempt];
end
writetable(h,'../data/tidy_data/lognormal.csv')

