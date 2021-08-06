import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statistics as stats
import math
def arreglardatos(datos):
     indices = datos.columns
     datoscorregidos = datos.copy()
     for column in indices:
         datoscorregidos[column] = datoscorregidos[column].str.replace("%","").astype(float)
     return datoscorregidos

def Rayleighdist(x,sigma):
    N= len(x)
    y = np.zeros(N)
    for k in range(N):
        if x[k]<0:
            y[k]= 0
        else:
            y[k]=((x[k])/sigma**2)*math.exp(-(x[k])**2/(2*sigma**2))
    return y

def gammadist(x,k,beta):
    N= len(x)
    y = np.zeros(N)
    for i in range(N):
        if x[i]<0:
            y[i]= 0
        else:
            y[i]=((beta**(k))*(x[i]**(k-1))/(math.gamma(k)))*math.exp(-beta*x[i])
    return y

#Lectura de datos
data = pd.read_csv('data/raw_data/M-0319-2017.csv',  header=0,sep = ',',encoding='ISO-8859-1', low_memory=False)
porcentajes = data[['Phase C-A V-Harmonic 3rd', 'Phase C-A V-Harmonic 5th']]
ca = [x for x in data.columns if 'Phase C-A V-Harmonic' in x]
bc = [x for x in data.columns if 'Phase B-C V-Harmonic' in x] 
ab = [x for x in data.columns if 'Phase A-B V-Harmonic' in x] 
porcentajes= data[ca + bc + ab]

#Generación de archivo para análisis con EM
datoscorregidos = arreglardatos(porcentajes)
datoscorregidos.to_csv('data/tidy_data/datos.csv', index = False)

#Generación de histogramas

for variable in datoscorregidos.columns:
    plt.figure()
    plt.hist(arreglardatos(porcentajes)[variable], bins = 'auto', density = True)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.savefig('figures/exploratory_figures/hist' + variable.replace(" ","").replace("-","_")+'.png')
    plt.close('all')   
#Exploración de parámetros estadísticos
pares = arreglardatos(porcentajes).mean()[1::2] 
impares = arreglardatos(porcentajes).mean()[::2]
print('Promedio pares', pares.mean())
print('Promedio impares', impares.mean())
paresca = datoscorregidos[ca].mean()[1::2] 
imparesca = datoscorregidos[ca].mean()[::2]
print('Promedio pares ca', paresca.mean())
print('Promedio impares ca', imparesca.mean())
paresbc = datoscorregidos[bc].mean()[1::2] 
imparesbc = datoscorregidos[bc].mean()[::2]
print('Promedio pares bc', paresbc.mean())
print('Promedio impares bc', imparesbc.mean())
paresab = datoscorregidos[ab].mean()[1::2] 
imparesab = datoscorregidos[ab].mean()[::2]
print('Promedio pares ab', paresab.mean())
print('Promedio impares ab', imparesab.mean())
estadisticos = pd.DataFrame(ca + bc + ab, columns = ['Variable'])
estadisticos['Promedio'] = arreglardatos(porcentajes).mean().reset_index(drop=True)
estadisticos['Desviación estándar'] = arreglardatos(porcentajes).std().reset_index(drop=True)
#Lectura de resultados del análsis
gamma = pd.read_csv('data/tidy_data/gamma.csv',  header=0,sep = ',',encoding='ISO-8859-1', low_memory=False)
ray = pd.read_csv('data/tidy_data/rayleigh.csv',  header=0,sep = ',',encoding='ISO-8859-1', low_memory=False)
#Generación de histogramas con gráficos de resultados
for i, variable in enumerate(datoscorregidos.columns[1::2]):
    try:
        plt.figure()
        plt.hist(datoscorregidos[variable], bins = 'auto', density = True)
        x0 = datoscorregidos[variable].min()
        x1 = datoscorregidos[variable].max()
        x = np.arange(0,x1*1.1,0.00001)
        y1 = ray['alfa1'][i]*Rayleighdist(x, ray['sigma1'][i])+ray['alfa2'][i]*Rayleighdist(x, ray['sigma2'][i])+ray['alfa3'][i]*Rayleighdist(x, ray['sigma3'][i])
        y2 = gamma['alfa1'][i]*gammadist(x,gamma['k1'][i],gamma['beta1'][i]) + gamma['alfa2'][i]*gammadist(x,gamma['k2'][i],gamma['beta2'][i])+gamma['alfa3'][i]*gammadist(x,gamma['k3'][i],gamma['beta3'][i])
        plt.plot(x,y1,label = 'Rayleigh')
        plt.plot(x,y2, label = 'Gamma')
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.legend()
        plt.savefig('figures/explanatory_figures/hist' + variable.replace(" ","").replace("-","_")+'.png')
        plt.close('all') 
    except:
        print('No se pudo analizar' + variable)
