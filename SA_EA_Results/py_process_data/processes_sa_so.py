# -*- coding: utf-8 -*-
"""
Created on Wed May  5 14:27:05 2021

@author: yl918888
"""

import os 
import numpy as np
import pandas as pd
#%%
root = r"C:\Users\yl918888\Desktop\SA_EA_SO_EXP\SO_Results_NF_10K"
listdir = os.listdir(root)
#%%

DE = ["method", "algo", "problem", "pop_size", "base_vector", "diff_vector-pop_size_ratio", "crossover_method", "crossover_prob", "beta_min", "beta_ext", "BestSol_mean", "BestSol_std"]
CMAES = ["method", "algo", "problem", "lambda", "mu-lambda_ratio", "sigma0", "alpha_mu", "rescale_sigma0", "BestSol_mean", "BestSol_std"]
Algo = ['DE'] # 'DE','CMA-ES'

for algo in Algo:
    print(algo)
    if algo == "DE":
        dfAlgo = pd.DataFrame(columns = DE)
    else:
        dfAlgo = pd.DataFrame(columns = CMAES)
    
    for dirs in listdir:
        if ('Problem' in dirs and algo in dirs):
            listFiles = os.listdir(os.path.join(root,os.path.join(dirs,'ModelEval')))
            dirsVal = dirs.split(' ')[2]
            print('  ',dirs)
            for files in listFiles:
                print('    ',files)
                params = [pVal for pVal in dirsVal.split('-') if 'ES' not in pVal]
                paramsTemp = [float(param) for param in files.split('-') if param not in 'metrics.csv']
                params.extend(paramsTemp)
                pathData = os.path.join(root,os.path.join(os.path.join(dirs,'ModelEval'),files))
                data = pd.read_csv(pathData)
                BestSolMean = data.mean()['BestSolution']
                BestSolStd = data.std()['BestSolution']
                params.append(BestSolMean)
                params.append(BestSolStd)
                if algo == "DE":
                    dfTemp = pd.DataFrame([params], columns=DE)
                else: #
                    dfTemp = pd.DataFrame([params], columns=CMAES)
                dfAlgo = pd.concat([dfAlgo,dfTemp])
            # end for all files
        #end if dir problem 
    #end all algo
    dfAlgo.to_csv(algo+"_All.csv", index=False)

                