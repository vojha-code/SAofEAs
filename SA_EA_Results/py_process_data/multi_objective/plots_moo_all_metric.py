
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 14:27:05 2021

@author: yl918888
"""

import os 
import numpy as np
import pandas as pd
#%%
root = r"C:\Users\yl918888\Desktop\Results_Collection_SA_EA\Collection_Results_SA_MOO"
listdir = os.listdir(root)
#%%

## "pop", "pSBX", "sbxDI"", "pm", "pmDI", "tournament"
NSGAIII = ["method", "algo", "problem","dim", "pop", "pSBX", "sbxDI", "pm", "pmDI", "tournament", "GD", "IGD", "HV"]
##"pop", "pSBX", "sbxDI", "pm", "pmDI", "mode", "neighbour"
MOEAD = ["method", "algo", "problem", "dim", "pop", "pSBX", "sbxDI", "pm", "pmDI", "mode", "neighbour",  "GD", "IGD", "HV"]

Algo = ['NSGAIII', 'MOEAD'] # 'DE','CMA-ES'

for algo in Algo:
    print(algo)
    if algo == "NSGAIII":
        dfAlgo = pd.DataFrame(columns = NSGAIII)
    else:
        dfAlgo = pd.DataFrame(columns = MOEAD)
    
    for dirs in listdir:
        if ('Problem' in dirs and algo in dirs):
            listFiles = os.listdir(os.path.join(root,os.path.join(dirs,'ModelEval')))
            dirsVal = dirs.split(' ')[2]
            print('  ',dirs)
            for files in listFiles:
                print('    ',files)
                params = [pVal for pVal in dirsVal.split('-') if 'ES' not in pVal]
                params.remove('Problem')
                paramsTemp = [float(param) for param in files.split('-') if param not in 'metrics.csv']
                params.extend(paramsTemp)
                pathData = os.path.join(root,os.path.join(os.path.join(dirs,'ModelEval'),files))
                data = pd.read_csv(pathData)
                Best_GD_Mean = data.mean()['GD']
                Best_IGD_Mean = data.std()['IGD']
                Best_HV_Mean = data.std()['HV']
                params.append(Best_GD_Mean)
                params.append(Best_IGD_Mean)
                params.append(Best_HV_Mean)
                if algo == "NSGAIII":
                    dfTemp = pd.DataFrame([params], columns=NSGAIII)
                else: #
                    dfTemp = pd.DataFrame([params], columns=MOEAD)
                dfAlgo = pd.concat([dfAlgo,dfTemp])
                
            # end for all files
        #end if dir problem 
    #end all algo
    dfAlgo.to_csv(algo+"_All.csv", index=False)









