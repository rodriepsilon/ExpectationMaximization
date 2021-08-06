classdef gammaMix
    properties
        NumComponents
        LogLikelihood
        ComponentProportions
        k
        beta
    end
    methods
        function obj = gammaMix(N, Likelihood ,alpha ,k, beta) 
            %se define un nuevo objeto de la clase
            obj.NumComponents = N;
            obj.LogLikelihood = Likelihood ;
            obj.ComponentProportions = alpha ;
            obj.k = k;
            obj.beta = beta;
        end
        function f = pdf(obj ,X) 
            %se calcula la mezcla de densidades en un vector de datos
            f = 1:length (X);
            for i = 1:length(X)
                f(i) = 0;
                for n = 1:obj.NumComponents 
                    %se calcula en cadacomponente y se multiplica por el peso
                    val = Gamma(X(i),obj.k(n),obj.beta(n));
                    f(i) = f(i) + obj.ComponentProportions (n)*val;
                end
            end
        end
        function r = random (obj ,n) 
            %se genera un vector "n" de valores aleatorios distribuidos 
            %de acuerdo con el ajuste
            priori = zeros (1, obj.NumComponents + 1);
            for k = 1:obj. NumComponents
                for j = 1:i 
  %se separa el intervalo [0 ,1] en segmentos que corresponden a 
  %la probabilidad de escoger un componente
                    priori(i) = priori(i)+obj.ComponentProportions(j);
                end
                r = zeros (1,n);
                for i =1:n
                    dr = rand; 
 %este valor aleatorio uniforme indica el componente de donde proviene el dato
                    actComp = 1;
                    while priori(actComp)<dr
                        actComp = actComp + 1;
                    end
                    r(i) = random('Gamma',obj.k( actComp ),1/obj.beta(actComp));
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
                K = 1:N;
                Beta = 1:N;
                for rep = 1:5 %cantidad de veces que se aplica el algoritmo
                    for i = 1:N 
            %los valores iniciales se escogen aleatoriamente en intervalos razonables
                        Alpha (i) =1/N;
                        K(i) =((Mean^2)/Variance)+(0.0 + 2.0* rand );
                        Beta (i) = (Mean/Variance)+(0.0 + 2.0*rand );
                    end
                    for t = 1:500 
                       %cantidad de iteraciones cada vez que se aplica algoritmo
                        AlphaNew = Alpha;
                        KNew = K;
                        BetaNew = Beta;
                        ExpectedValueMatrix = zeros (N,M); %matriz de valores esperados de los datos invisibles
                    for i = 1:N
                        for m = 1:M
                            ExpectedValueMatrix (i,m) = ExpectedValueYGamma(i,m,X,Alpha,K,Beta);
                        end
                    end
                    for i = 1:N %se calculan los nuevos valores de " Alpha " y "Mu"
                        sum1 = 0;
                        sum2 = 0;
                        sum3 = 0;
                        for m = 1:M
                            sum1 = sum1 + ExpectedValueMatrix (i,m);
                            sum2 = sum2 + ExpectedValueMatrix (i,m) * X(m);
                        end
                        
                        for m = 1:M
                            sum3 = sum3 + ExpectedValueMatrix(i,m)*(log(sum2/sum1)-log(X(m)));
                        end
                        AlphaNew (i) = sum1 /M;
                        KNew (i) = sum1 /(2*sum3);
                        if isnan( AlphaNew(i)) 
                            %en caso deerror , se mantiene el valoranterior
                            AlphaNew (i) = Alpha(i);
                        end
                        if isnan ( KNew (i)) 
                            %en caso de error ,se mantiene el valor anterior
                            KNew (i) = K(i);
                        end
                    end
                    for i = 1:N %se calculan los nuevos valores de "beta"
                        sum1 = 0;
                        sum2 = 0;
                        for m = 1:M
                            sum1 = sum1 + ExpectedValueMatrix (i,m)*X(m);
                            sum2 = sum2 + ExpectedValueMatrix (i,m);
                        end
                        BetaNew (i) = (sum2*KNew(i))/sum1;
                        if BetaNew (i)==0 || isnan ( BetaNew (i)) 
                            %en caso de error , se mantiene elvalor anterior
                            BetaNew (i) = Beta (i);
                        end
                    end
                    if EliminateImpulses %se verifica si esnecesario eliminar impulsos
                        [maxKNew , maxKNewIndex] = max(KNew ); 
                        %se encuentra el mayor "Sigma " y su lugar en el vector
                        for i =1:N
                            if maxKNew > 1000* KNew (i)
                %si el " Sigma " es mucho menor ,se cambia por uno razonable
                                KNew (i) = KNew (maxKNewIndex );
                                BetaNew (i) = BetaNew(maxKNewIndex);
                            end
                        end
                    end
                    if t == 1
                        Qprev = ExpectedValueQGamma(X, Alpha ,K,Beta, AlphaNew ,KNew , BetaNew );
                    else
                        Qnext = ExpectedValueQGamma(X, Alpha ,K ,Beta, AlphaNew ,KNew , BetaNew );
                        if Qnext *1.01 < Qprev 
                            %si no hubo mejoras significativas en el valor
                            break
                        end
                        Qprev = Qnext ;
                    end
                    Alpha = AlphaNew ;
                    K = KNew ;
                    Beta = BetaNew ;
                end
                Likelihood = LogLikelihoodGamma(X,Alpha ,K, Beta);
                if rep ==1
                    LikelihoodFinal = Likelihood ;
                    AlphaFinal = Alpha ;
                    KFinal = K;
                    BetaFinal = Beta;
                elseif Likelihood > LikelihoodFinal 
       %si los nuevos datos son mejores , estos se guardan como valores finales
                    if not( isnan ( Likelihood ) || isnan (LikelihoodFinal ))
                        LikelihoodFinal = Likelihood ;
                        AlphaFinal = Alpha;
                        KFinal = K;
                        BetaFinal = Beta;
                    end
                end
            end
            obj = gammaMix(N, LikelihoodFinal , AlphaFinal, KFinal,BetaFinal);
        end
    end
end