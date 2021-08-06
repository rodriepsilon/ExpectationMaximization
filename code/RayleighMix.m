classdef RayleighMix
    properties
        NumComponents
        LogLikelihood
        ComponentProportions
        sigma
    end
    methods
        function obj = RayleighMix(N, Likelihood ,alpha, sigma)
            %se define un nuevo objeto de la clase
            obj.NumComponents = N;
            obj.LogLikelihood = Likelihood ;
            obj.ComponentProportions = alpha ;
            obj.sigma = sigma ;
        end
        function f = pdf(obj ,X) %se calcula la mezcla de densidades en un vector de datos
            f = 1:length (X);
            for k = 1:length(X)
                f(k) = 0;
                for n = 1:obj.NumComponents 
                    %se calcula en cadacomponente y se multiplica por el peso
                    val = Rayleigh(X(k),obj.sigma(n));
                    f(k) = f(k) + obj.ComponentProportions(n)*val;
                end
            end
        end
        function r = random(obj ,n)
            %se genera un vector "n" de valores aleatorios 
            %distribuidos de acuerdo con el ajuste
            priori = zeros (1, obj.NumComponents + 1);
            for k = 1:obj. NumComponents
                for j = 1:k
               %se separa el intervalo [0 ,1] en segmentos que 
               %corresponden a la probabilidad de escoger un componente
                    priori(k) = priori(k)+obj.ComponentProportions(j);
                end
                r = zeros(1,n);
                for k =1:n
                    dr = rand; 
                    %este valor aleatorio uniforme indica el componente de donde proviene el dato
                    actComp = 1;
                    while priori(actComp)<dr
                        actComp = actComp + 1;
                    end
                    r(k) = random('Rayleigh',obj.sigma( actComp ));
                end
            end
        end
    end
        methods(Static)
            function obj = fit(X,N, EliminateImpulses ) 
                %se genera el ajuste a los datos "X" utilizando "N" componentes
                if nargin ==2
                    EliminateImpulses = true ; %por defecto , se eliminan los impulsos
                end
                M = length (X);
                Mean = mean (X);
                Variance = var(X);
                Alpha = 1:N;
                %Mu = 1:N;
                Sigma = 1:N;
                for rep = 1:5 %cantidad de veces que se aplica el algoritmo
                    for k = 1:N 
                        %los valores iniciales se escogen aleatoriamente en intervalos razonables
                        %Valores iniciales
                        Alpha (k) =1/N;                      
                        Sigma (k) = sqrt (Mean*sqrt(2/pi))* (0.0 + 2*rand );
                    end
                    for t = 1:50 %cantidad de iteraciones cada vez que se aplica algoritmo
                        AlphaNew = Alpha ;
                        SigmaNew = Sigma ;
                        ExpectedValueMatrix = zeros (N,M); 
                        %matriz de valores esperados de los datos invisibles
                    for i = 1:N
                        for m = 1:M
                            ExpectedValueMatrix(i,m) = ExpectedValueYRayleigh(i,m,X,Alpha, Sigma );
                        end
                    end
                    for i = 1:N %se calculan los nuevos valores de " Alpha " 
                        sum1 = 0;
                        for m = 1:M
                            sum1 = sum1 + ExpectedValueMatrix (i,m);                      
                        end
                        AlphaNew (i) = sum1 /M;
                        if isnan ( AlphaNew (i)) 
                            %en caso deerror , se mantiene el valoranterior
                            AlphaNew (i) = Alpha (i);
                        end
                    end
                    for i = 1:N %se calculan los nuevos valores de " Sigma "
                        sum = 0;
                        for m = 1:M
                            sum = sum + ExpectedValueMatrix (i,m)*X(m)^2;
                        end
                        SigmaNew (i) = sqrt ( sum/ ( AlphaNew (i)* M*2) );
                        if SigmaNew (i)==0 || isnan ( SigmaNew (i)) %en caso de error , se mantiene elvalor anterior
                            SigmaNew (i) = Sigma (i);
                        end
                    end
                    if EliminateImpulses
                        %se verifica si esnecesario eliminar impulsos
                        [maxSigmaNew , maxSigmaNewIndex] = max(SigmaNew ); 
                        %se encuentra el mayor "Sigma " y su lugar en el vector
                        for i =1:N
                            if maxSigmaNew > 1000* SigmaNew (i)
                      %si el " Sigma " es mucho menor ,se cambia por uno razonable
                                SigmaNew (i) = SigmaNew (maxSigmaNewIndex );
                            end
                        end
                    end
                    if t == 1
                        Qprev = ExpectedValueQRayleigh(X, Alpha ,Sigma , AlphaNew, SigmaNew );
                    else
                        Qnext = ExpectedValueQRayleigh(X, Alpha,Sigma , AlphaNew , SigmaNew );
                        if Qnext *1.01 < Qprev 
                            %si no hubo mejoras significativas en el valor
                            break
                        end
                        Qprev = Qnext ;
                    end
                    Alpha = AlphaNew ;
                    Sigma = SigmaNew ;
                end
                Likelihood = LogLikelihoodRayleigh(X,Alpha, Sigma );
                if rep ==1
                    LikelihoodFinal = Likelihood ;
                    AlphaFinal = Alpha ;
                    SigmaFinal = Sigma ;
                elseif Likelihood > LikelihoodFinal 
                    %si los nuevos datos son mejores , estos se guardan como valores finales
                    if not( isnan ( Likelihood ) || isnan (LikelihoodFinal ))
                        LikelihoodFinal = Likelihood ;
                        AlphaFinal = Alpha ;
                        %MuFinal = Mu;
                        SigmaFinal = Sigma ;
                    end
                end
            end
            obj = RayleighMix (N, LikelihoodFinal , AlphaFinal , SigmaFinal );
        end
    end
end