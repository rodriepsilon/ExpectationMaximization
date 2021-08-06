from sklearn import mixture
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('data/tidy_data/datos.csv',  header=0,sep = ',',encoding='ISO-8859-1', low_memory=False)
gmm = mixture.GaussianMixture(n_components=3, max_iter=1000, covariance_type='full').fit(data['Phase B-C V-Harmonic 3rd'].values.reshape(-1,1))
'''
Código para comparar resultados de log-verosimilitud de SkLearn con el modelo de MatLab
'''
print('Loglikelihood')
print(gmm.score_samples(data['Phase B-C V-Harmonic 3rd'].values.reshape(-1,1)).sum())
