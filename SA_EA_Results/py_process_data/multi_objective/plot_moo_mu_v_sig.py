# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 11:51:26 2021

@author: yl918888
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 16:20:46 2021


@author: yl918888
"""
import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plots = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\py_process_data\plots'
folder = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_MOO'
#%%

##pop	pSBX	sbxDI	pm	pmDI	mode	neighbour
##pop	pSBX	sbxDI	pm	pmDI	tournament

MOEAD = [r"$\lambda$", r"$P[\mathrm{X}]$", r"$\mathrm{X}_{\mathrm{DI}}$", r"$P[\mathrm{PM}]$", r"$\mathrm{PM}_{\mathrm{DI}}$",r"$Mode$", r"$\epsilon_N$"]
NSGAIII = [r"$\lambda$", r"$P[\mathrm{X}]$", r"$\mathrm{X}_{\mathrm{DI}}$", r"$P[\mathrm{PM}]$", r"$\mathrm{PM}_{\mathrm{DI}}$", r"$K$"]

#%%
markerList = ["o", "v", "*", "X", "D", "+", "3", "s"]

AlgoList = ['moead','nsgaiii']

for algo in AlgoList:
    root = os.path.join(folder,algo)
    listdir = os.listdir(root)
   
    fig, axs = plt.subplots(3, 3, figsize=(7, 8), sharey=True, sharex=True)
    
    checkTest = 'var'
    i = 0
    j = 0
    for file in listdir:
        if checkTest in file:
            print(file)
            pathData = os.path.join(root,file)
            data = pd.read_csv(pathData)
            columns = data.columns.tolist()
            
            columnsMu = [colName for colName in columns if '.1' not in colName]
            columnsSi = [colName for colName in columns if '.1' in colName]
            dataMu = data[columnsMu]
            dataSi = data[columnsSi]
            dataMuT = dataMu.T       
            dataMuTnorm=(dataMuT-dataMuT.min())/(dataMuT.max()-dataMuT.min())
            dataSiT = dataSi.T       
            dataSiTnorm=(dataSiT-dataSiT.min())/(dataSiT.max()-dataSiT.min())
            
            dataMunorm = dataMuTnorm.T
            dataSinorm = dataSiTnorm.T
            
            cmap = plt.cm.get_cmap('Accent') 
            #if '_GD' in file:
            #    cmap = plt.cm.get_cmap('tab10') 
            #elif '_GD' in file:
            #    cmap = plt.cm.get_cmap('Accent') 
            #else:
            #    cmap = plt.cm.get_cmap('Accent') 
            
            #colors = [cmap(each) for each in np.linspace(0, 1, len(columnsMu))]
            colors = [cmap(each) for each in np.linspace(0, 1, 7)]

            for paramIndex in range(len(columnsMu)):
                paramX = columnsMu[paramIndex]
                paramY = columnsSi[paramIndex]
                if 'MOEAD' in file:
                    paramSTR = MOEAD[paramIndex]
                    plotAlgo = "MOEAD"
                else:
                    paramSTR = NSGAIII[paramIndex]
                    plotAlgo = "NSGAIII"
                     
                x = dataMunorm[paramX]
                y = dataSinorm[paramY]
                xmean = np.mean(dataMunorm[paramX].dropna().tolist())
                ymean = np.mean(dataSinorm[paramY].dropna().tolist())
            
                xstd = np.mean(dataMunorm[paramX].dropna().tolist())
                ystd = np.mean(dataSinorm[paramY].dropna().tolist())
                from matplotlib.patches import Ellipse
                ellipse = Ellipse(xy=(xmean, ymean), width=xstd, height=ystd, edgecolor=colors[paramIndex], fc='None', lw=1)
                axs[i,j].add_patch(ellipse)
                #axs[i,j].scatter(x, y, alpha=0.8, marker=markerList[paramIndex], label=paramSTR, color=colors[paramIndex])
                axs[i,j].scatter(xmean, ymean, alpha=1, marker=markerList[paramIndex], label=paramSTR, color=colors[paramIndex])
                #axs[i,j].plot((0,1), (1,1), alpha=0.1, color='gray', lw=0.5)
                #axs[i,j].plot((1,1), (1,0), alpha=0.1, color='gray', lw=0.5)
                
                axs[i,j].set_ylim(-0.05,1.05)
                axs[i,j].set_xlim(-0.05,1.05)
            
            if i == 0 and j == 0:
                axs[i,j].set_title('GD')
            if i == 0 and j == 1:
                axs[i,j].set_title('IGD')
            if i == 0 and j == 2:
                axs[i,j].set_title('HV')
    
            if i == 0 and j == 0:
                axs[i,j].set_ylabel("Morris LHS  $\sigma$")          
            if i == 0:
                axs[i,j].set_xlabel("Morris LHS  $\mu$") 
            
            if i == 1 and j == 0:
                axs[i,j].set_ylabel("Morris  $\sigma$")          
            if i == 1:
                axs[i,j].set_xlabel("Morris  $\mu$") 
                
            if i == 2:
                axs[i,j].set_xlabel("Sobol $S_i$")         
                
            if i == 2 and j == 0: 
                axs[i,j].set_ylabel("Sobol  $ST_i$")  
                if 'MOEAD' in file:
                    axs[i,j].legend(bbox_to_anchor=(0.5, 0), labelspacing=0.2, loc="lower center",  bbox_transform=fig.transFigure, ncol=7)
                else:
                    axs[i,j].legend(bbox_to_anchor=(0.5, 0), labelspacing=1, loc="lower center",  bbox_transform=fig.transFigure, ncol=6)

            # GO NEXT PLOT            
            j = j + 1
            if j%3 == 0:
                j = 0
                i = i + 1
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots, "Plot_Scatter_"+plotAlgo+".pdf"),bbox_inches='tight')
    plt.savefig(os.path.join(plots, "Plot_Scatter_"+plotAlgo+".png"),bbox_inches='tight')
    plt.show()

#%% Bar plot

def sortValue(xAll, yAll, paramSTR):
    x =  np.zeros(len(xAll[0]))
    y =  np.zeros(len(yAll[0]))
    #for k in range(len(xAll)):
    for k in range(3):        
        x = x + xAll[k]
        y = y + yAll[k]
    x =  (x - np.min(x)) / (np.max(x) - np.min(x))
    y =  (y - np.min(y)) / (np.max(y) - np.min(y))            
    # Sum            
    xpy = x+y
    indexSort = np.argsort(xpy)
        
    paramSort = [paramSTR[idx] for idx in indexSort]
    xSort = [x[idx] for idx in indexSort]
    ySort = [y[idx] for idx in indexSort]
    return paramSort, xSort, ySort

AlgoList = ['moead','nsgaiii']

MetricList = ['_GD_', '_IGD_', '_HV_']
SampleList = ['MorisLHS_', 'Moris_', 'SOBAL_']

for algo in AlgoList:
    root = os.path.join(folder,algo)
    listdir = os.listdir(root)
    fig, axs = plt.subplots(4, 4, figsize=(12, 6), sharex=False, sharey=True)
    
    
    xAll_GD = []
    xAll_IGD = []
    xAll_HV = []
    
    yAll_GD = []
    yAll_IGD = []
    yAll_HV = []

    xAll_L = []
    xAll_M = []
    xAll_S = []
    
    yAll_L = []
    yAll_M = []
    yAll_S = []    
    
    checkTest = 'var'
    i = 0
    j = 0
    for file in listdir:
        if checkTest in file:
            pathData = os.path.join(root,file)
            data = pd.read_csv(pathData)
            columns = data.columns.tolist()
            columnsMu = [colName for colName in columns if '.1' not in colName]
            columnsSi = [colName for colName in columns if '.1' in colName]
            dataMu = data[columnsMu]
            dataSi = data[columnsSi]
            
            dataMuT = dataMu.T       
            dataMuTnorm=(dataMuT-dataMuT.min())/(dataMuT.max()-dataMuT.min())
            dataSiT = dataSi.T       
            dataSiTnorm=(dataSiT-dataSiT.min())/(dataSiT.max()-dataSiT.min())
            
            dataMunorm = dataMuTnorm.T
            dataSinorm = dataSiTnorm.T
            
            dataMuNormSum = dataMunorm.sum()
            dataSiNormSum = dataSinorm.sum()
                
            if 'MOEAD' in file:
                paramSTR = MOEAD
                plotAlgo = "MOEAD"
            else:
                paramSTR = NSGAIII
                plotAlgo = "NSGAIII"
    
            if j < 2:             
                xLabel = "$\mu$"
                yLabel = "$\sigma$"
            else:
                xLabel = "$S_i$"
                yLabel = "$ST_i$"
                    
                
            #colors = [cmap(each) for each in np.linspace(0, 1, len(columnsMu))]      
            colors = ['black', 'lightgray']
            width = 0.35
            x = dataMuNormSum.tolist()
            y = dataSiNormSum.tolist()
            x =  np.asarray(x)
            x =  (x - np.min(x)) / (np.max(x) - np.min(x))
            y =  np.asarray(y)
            y =  (y - np.min(y)) / (np.max(y) - np.min(y))
            
            if '_GD_' in file:
                xAll_GD.append(x)    
                yAll_GD.append(y)
                #print('GD >>>',i,j)
            
            if '_IGD_' in file:
                xAll_IGD.append(x)    
                yAll_IGD.append(y)    
                
            if '_HV_' in file:
                xAll_HV.append(x)    
                yAll_HV.append(y)    
                
            if 'MorisLHS_' in file:
                #print('     ', file)
                xAll_L.append(x)    
                yAll_L.append(y)    
                
            if 'Moris_' in file:
                xAll_M.append(x)    
                yAll_M.append(y)    
                
            if 'SOBAL_' in file:
                xAll_S.append(x)    
                yAll_S.append(y)    
            
            xpy = x+y
            indexSort = np.argsort(xpy)
            
            paramSort = [paramSTR[idx] for idx in indexSort]
            xSort = [x[idx] for idx in indexSort]
            ySort = [y[idx] for idx in indexSort]
        
            #axs[i,j].bar(paramSTR, x, width, label =xLabel, color=colors[0], alpha=0.7)
            #axs[i,j].bar(paramSTR, y, width, bottom = x, label=yLabel, color=colors[1],alpha=0.7)
            
            axs[i,j].bar(paramSort, xSort, width, label =xLabel, color="w", alpha=1, edgecolor=colors[0], linewidth=1)
            axs[i,j].bar(paramSort, ySort, width, bottom = xSort, label=yLabel, color=colors[1],alpha=1, edgecolor=colors[0], linewidth=1)
            print(i,j,algo,file)
            
            #if j == 0:
                #print(j, '>>>>>>>', file)
            
            if i == 0:                
                print(i, '>>>>>>>', file)
            if j == 2 and i == 2:
                axs[i,j].legend(ncol=2, loc='upper left')
            else:
                axs[i,j].legend(ncol=2)
    
            if i == 0 and j == 0:
                axs[i,j].set_ylabel('GD')
            if i == 1 and j == 0:
                axs[i,j].set_ylabel('IGD')
            if i == 2 and j == 0:
                axs[i,j].set_ylabel('HV')
                
            if i == 0 and j == 0:
                axs[i,j].set_title('Morris LHS')
            if i == 0 and j == 1:
                axs[i,j].set_title('Morris')
            if i == 0 and j == 2:
                axs[i,j].set_title('Sobol')
            
            # GO NEXT PLOT            
            i = i + 1
            if i%3 == 0:
                i = 0
                j = j + 1
 
    axs[3,0].set_ylabel('GD + IGD + HV')
    axs[0,3].set_title('LHS + Morris + Sobol')
    
    xLabel = "direct"
    yLabel = "interaction"
    
    for k in range(3):
        if k == 0:
            paramSort_Met, xSort_Met, ySort_Met = sortValue(xAll_GD, yAll_GD, paramSTR)
            paramSort_Sam, xSort_Sam, ySort_Sam = sortValue(xAll_L, yAll_L, paramSTR)
        if k == 1:
            paramSort_Met, xSort_Met, ySort_Met = sortValue(xAll_IGD, yAll_IGD, paramSTR)
            paramSort_Sam, xSort_Sam, ySort_Sam = sortValue(xAll_M, yAll_M, paramSTR)
        if k == 2:
            paramSort_Met, xSort_Met, ySort_Met = sortValue(xAll_HV, yAll_HV, paramSTR)
            paramSort_Sam, xSort_Sam, ySort_Sam = sortValue(xAll_S, yAll_S, paramSTR)          
            
        
        axs[k,3].bar(paramSort_Met, xSort_Met, width, label =xLabel, color="w", alpha=1, edgecolor=colors[0], linewidth=1)
        axs[k,3].bar(paramSort_Met, ySort_Met, width, bottom = xSort_Met, label=yLabel, color=colors[1],alpha=1, edgecolor=colors[0], linewidth=1)
        axs[k,3].legend(ncol=1, loc='upper left')
        print(k,3,algo,file)
        
        axs[3,k].bar(paramSort_Sam, xSort_Sam, width, label =xLabel, color="w", alpha=1, edgecolor=colors[0], linewidth=1)
        axs[3,k].bar(paramSort_Sam, ySort_Sam, width, bottom = xSort_Sam, label=yLabel, color=colors[1],alpha=1, edgecolor=colors[0], linewidth=1)
        axs[3,k].legend(ncol=1, loc='upper left')
        print(k,3,algo, file)
        
    #axs[3,3].plot([], [], ' ')
    axs[3,3].set_visible(False)
    
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots, "Plot_Bar_"+plotAlgo+".pdf"),bbox_inches='tight')
    plt.savefig(os.path.join(plots, "Plot_Bar_"+plotAlgo+".png"),bbox_inches='tight')
    plt.show()

#%% Line plot

#MOEAD = [r"$\lambda$", r"$P[\mathrm{X}]$", r"$\mathrm{X}_{\mathrm{DI}}$", r"$P[\mathrm{PM}]$", r"$\mathrm{PM}_{\mathrm{DI}}$",r"$Mode$", r"$\epsilon_N$"]
#NSGAIII = [r"$\lambda$", r"$P[\mathrm{X}]$", r"$\mathrm{X}_{\mathrm{DI}}$", r"$P[\mathrm{PM}]$", r"$\mathrm{PM}_{\mathrm{DI}}$", r"$K$"]


AlgoList = ['MOEAD', 'NSGAIII']
for algo in AlgoList:
    listdir = os.listdir(folder)
    
    checkTest = algo+'_All'
    i = 0
    j = 0
    for file in listdir:
        if checkTest in file:
            print(file)
            pathData = os.path.join(folder,file)
            data = pd.read_csv(pathData)    
            columns = data.columns.tolist()
            fig, axs = plt.subplots(3, len(columns)-7, figsize=(2*(len(columns)-4), (len(columns)-4)/2+3), sharey =True)
            #fig, axs = plt.subplots(1, len(columns)-4, sharey =True)
            
            if 'MOEAD' in file:
                cmap = plt.cm.get_cmap('nipy_spectral') 
                paramSTR = MOEAD
            else:
                cmap = plt.cm.get_cmap('nipy_spectral') 
                paramSTR = NSGAIII
    
            colors = [cmap(each) for each in np.linspace(0, 1, len(paramSTR))]                  
            Metric = ['GD', 'IGD', 'HV']
            for i in range(len(Metric)):
                j = 0
                for paramIndex in range(4,len(columns)-3):
                    panamLabel = paramSTR[paramIndex-4]
                    colorsVal = colors[paramIndex-4]
                    print(paramIndex,'  ',columns[paramIndex], '=', panamLabel)
                    paramX = columns[paramIndex]
   
                    dataXY = data[[paramX, Metric[i]]]
                    dataXYGroup = dataXY.groupby(paramX, as_index=False).mean()
            
                    x = dataXYGroup[paramX].tolist()
                    #y = dataXYGroup['BestSol_mean'].tolist()
                    bins = np.linspace(min(x), max(x), num=20).tolist()
                    groupsMean = dataXYGroup.groupby(pd.cut(dataXYGroup[paramX], bins)).mean()
                    groupsStd = dataXYGroup.groupby(pd.cut(dataXYGroup[paramX], bins)).std()
                    x = groupsMean[paramX].tolist()
                    #x = [k for k in range(len(groupsMean[paramX].tolist()))]
                    y = groupsMean[Metric[i]].tolist()
                    yStd = groupsStd[Metric[i]].tolist()
                    
                    #plt.scatter(x,y)
                    yVal =  np.asarray(y)
                    #if i  == 2: #or i == 1 :# HV and IGD revirsed
                    #   yNorm = (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))
                    #else:
                    yNorm = (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))
                        
                    axs[i,j].plot(x,yNorm, label=panamLabel, color=colorsVal)
                    axs[i,j].legend(fontsize=16)
        
                    j = j + 1
                
                axs[i,0].set_ylabel(Metric[i]+'  Score', fontsize=16)
                
            #plt.tight_layout()
            #plt.savefig(os.path.join(plots, file+".pdf"),bbox_inches='tight')
            #plt.savefig(os.path.join(plots, file+".png"),bbox_inches='tight')
            plt.show()

#%% Line plot Methods specific (trail)

#MOEAD = [r"$\lambda$", r"$P[\mathrm{X}]$", r"$\mathrm{X}_{\mathrm{DI}}$", r"$P[\mathrm{PM}]$", r"$\mathrm{PM}_{\mathrm{DI}}$",r"$Mode$", r"$\epsilon_N$"]
#NSGAIII = [r"$\lambda$", r"$P[\mathrm{X}]$", r"$\mathrm{X}_{\mathrm{DI}}$", r"$P[\mathrm{PM}]$", r"$\mathrm{PM}_{\mathrm{DI}}$", r"$K$"]
from matplotlib import colors as mcolors
from scipy.ndimage import gaussian_filter1d
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

AlgoList = ['MOEAD', 'NSGAIII']
for algo in AlgoList:
    listdir = os.listdir(folder)
    
    checkTest = algo+'_All'
    i = 0
    j = 0
    for file in listdir:
        if checkTest in file:
            print(file)
            pathData = os.path.join(folder,file)
            data = pd.read_csv(pathData)    
            columns = data.columns.tolist()
            fig, axs = plt.subplots(3, len(columns)-7, figsize=(2*(len(columns)-4), (len(columns)-4)/2+3), sharey =True)
            #fig, axs = plt.subplots(1, len(columns)-4, sharey =True)
            
            if 'MOEAD' in file:
                #cmap = plt.cm.get_cmap('nipy_spectral') 
                paramSTR = MOEAD
                plotName =  "MOEAD_All"
            else:
                #cmap = plt.cm.get_cmap('nipy_spectral') 
                paramSTR = NSGAIII
                plotName =  "NSGAIII_All"
    
            cmap = plt.cm.get_cmap('nipy_spectral') 
            #colors = [cmap(each) for each in np.linspace(0, 1, len(paramSTR))]                  
            colors = [cmap(each) for each in np.linspace(0, 1, 7)]                  
            Metric = ['GD', 'IGD', 'HV']
            #Metric = ['HV']
            for i in range(len(Metric)):
                j = 0
                print(Metric[i])
                for paramIndex in range(4,len(columns)-3):
                    panamLabel = paramSTR[paramIndex-4]
                    colorsVal = colors[paramIndex-4]
                    print(paramIndex,'  ',columns[paramIndex], '=', panamLabel, end='')
                    paramX = columns[paramIndex]
                    
                    dataXY = data[(data[Metric[i]] != 0)][[paramX, Metric[i]]]
                    dataXYGroup = dataXY.groupby(paramX, as_index=False).mean()
            
                    x = dataXYGroup[paramX].tolist()
                    #y = dataXYGroup['BestSol_mean'].tolist()
                    bins = np.linspace(min(x), max(x), num=20).tolist()
                    groupsMean = dataXYGroup.groupby(pd.cut(dataXYGroup[paramX], bins)).mean()
                    groupsMean = groupsMean.dropna()
                    groupsStd = dataXYGroup.groupby(pd.cut(dataXYGroup[paramX], bins)).std()
                    x = groupsMean[paramX].tolist()
                    #x = [k for k in range(len(groupsMean[paramX].tolist()))]
                    y = groupsMean[Metric[i]].tolist()
                    yStd = groupsStd[Metric[i]].tolist()
                    
                    #plt.scatter(x,y)
                    yVal =  np.asarray(y)
                    #if i  == 2: #or i == 1 :# HV and IGD revirsed
                    #   yNorm = (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))
                    #else:
                    yNorm = (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))
                    yNorm_smoothed = gaussian_filter1d(yNorm, sigma=1)
                    #yNorm = yVal
                        
                    #axs[i,j].plot(x,yNorm, label=panamLabel, color=colorsVal)
                    #axs[i,j].plot(x,yNorm,  color=colors[0], lw=2)
                    axs[i,j].plot(x,yNorm_smoothed,  color=colors[0], lw=2)
                    #axs[i,j].legend(fontsize=16)
                    
                    #colorsVal =  [colors['k'],colors['b'],colors['r']]
                    print('{0:.2f}'.format(min(y)), '{0:.2f}'.format(max(y)),)
                    Methods = ['MorrisLHS','Morris','SOBOL']#,
                    MethodLebel = ['l','m','s']
                    for k in range(len(Methods)):
                        dataXY = data[(data[columns[0]] == Methods[k]) & (data[Metric[i]] != 0)][[paramX, Metric[i]]]
                    
                        #dataXY = data[[paramX, Metric[i]]]
                        dataXYGroup = dataXY.groupby(paramX, as_index=False).mean()
                
                        x = dataXYGroup[paramX].tolist()
                        #y = dataXYGroup['BestSol_mean'].tolist()
                        bins = np.linspace(min(x), max(x), num=20).tolist()
                        groupsMean = dataXYGroup.groupby(pd.cut(dataXYGroup[paramX], bins)).mean()
                        groupsMean = groupsMean.dropna()
                        groupsStd = dataXYGroup.groupby(pd.cut(dataXYGroup[paramX], bins)).std()
                        x = groupsMean[paramX].tolist()
                        #x = [k for k in range(len(groupsMean[paramX].tolist()))]
                        y = groupsMean[Metric[i]].tolist()
                        yStd = groupsStd[Metric[i]].tolist()
                        
                        #plt.scatter(x,y)
                        yVal =  np.asarray(y)
                        #if i  == 2: #or i == 1 :# HV and IGD revirsed
                        #   yNorm = (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))
                        #else:
                        yNorm = (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))
                        yNorm_smoothed = gaussian_filter1d(yNorm, sigma=1)
                        #yNorm = yVal
                            
                        #axs[i,j].plot(x,yNorm, label=panamLabel+MethodLebel[k], color=colors[k+1])
                        #axs[i,j].legend(fontsize=16)
                        
                        #axs[i,j].plot(x,yNorm, alpha=0.5, color=colors[k+1])
                        axs[i,j].plot(x,yNorm_smoothed, alpha=0.5, color=colors[k+1])
                        #axs[i,j].plot(x,yNorm, color=colors[k], linestyle="-",marker="o")
                        if i == 2:
                            axs[i,j].set_xlabel(panamLabel, fontsize=16)
                    j = j + 1
                    
                
                axs[i,0].set_ylabel(Metric[i]+'  Score', fontsize=16)
                
            plt.tight_layout()
            plt.savefig(os.path.join(plots, plotName+".pdf"),bbox_inches='tight')
            plt.savefig(os.path.join(plots, plotName+".png"),bbox_inches='tight')
            plt.show()
##Experiment with smothing            
#from scipy.ndimage import gaussian_filter1d
#y_smoothed = gaussian_filter1d(yNorm, sigma=1)
#plt.plot(yNorm)
#plt.plot(y_smoothed)
#%% Cluster plot 

clusterMOEADk = [[2,3,2],
                 [3,2,4],
                 [2,3,4]]

clusterNSGAk = [[3,2,2],
                [3,2,3],
                [3,3,3]]

for algo in AlgoList:
    root = os.path.join(folder,algo)
    listdir = os.listdir(root)
    
    #plt.style.use('seaborn-whitegrid')
    fig, axs = plt.subplots(3, 3, figsize=(8, 8))#, sharey=True, sharex=True)
    
    checkTest = 'var'
    i = 0
    j = 0
    for file in listdir:
        if checkTest in file:
            print(file)
            pathData = os.path.join(root,file)
            data = pd.read_csv(pathData)
            columns = data.columns.tolist()
            columnsMu = [colName for colName in columns if '.1' not in colName]
            columnsSi = [colName for colName in columns if '.1' in colName]
            dataMu = data[columnsMu]
            dataSi = data[columnsSi]
            
            dataMuT = dataMu.T       
            dataMuTnorm=(dataMuT-dataMuT.min())/(dataMuT.max()-dataMuT.min())
            dataSiT = dataSi.T       
            dataSiTnorm=(dataSiT-dataSiT.min())/(dataSiT.max()-dataSiT.min())
            
            dataMuTnorm = dataMuTnorm.fillna(0)
            dataSiTnorm = dataSiTnorm.fillna(0)
            
            dataMunorm = dataMuTnorm.T
            dataSinorm = dataSiTnorm.T
            
            dataClusterFunc  = pd.concat([dataMunorm, dataSinorm], axis=1)
            dataClusterParams = pd.concat([dataMuTnorm, dataSiTnorm], axis=0)
            
            dataCluster = dataClusterFunc
            if 'MOEAD' in file:
                clusterK = clusterMOEADk
            else:
                clusterK = clusterNSGAk
                
            from sklearn.cluster import KMeans
            from sklearn.metrics import silhouette_samples, silhouette_score
            Sum_of_squared_distances = []
            silhouette_avg_values = []
            Xmat = dataCluster.to_numpy()
            K = range(2,len(dataCluster))
            for kIDX in K:
                km = KMeans(n_clusters=kIDX)
                km = km.fit(Xmat)
                #Sum_of_squared_distances.append(km.inertia_)
                cluster_labels = km.fit_predict(Xmat)
                silhouette_avg_values.append(silhouette_score(Xmat, cluster_labels))
            
            Sum_of_squared_distances = silhouette_avg_values
            np.argmax(silhouette_avg_values)
            #clusters_no = clusterK[i][j]
            clusters_no = 2+np.argmax(silhouette_avg_values)
            print(i,j,clusters_no) 

            kmeans = KMeans(n_clusters=clusters_no).fit(dataCluster)
            centroids = kmeans.cluster_centers_
            labels = kmeans.labels_
            
            from sklearn.decomposition import PCA
            pca = PCA(n_components=2)
            pca.fit(dataCluster.T)
            dataPCA =  pca.components_.T
            
            x =  dataPCA[:,0]
            y =  dataPCA[:,1]
                    
            if '_GD' in file:
                cmap = plt.cm.get_cmap('tab10') 
            elif '_IGD' in file:
                cmap = plt.cm.get_cmap('nipy_spectral') 
            else:
                cmap = plt.cm.get_cmap('Set1') 
                
            if 'MOEAD' in file:
                paramSTR = MOEAD
                plotAlgo = "MOEAD"
            else:
                paramSTR = NSGAIII
                plotAlgo = "NSGAIII"
                
            colors = [cmap(each) for each in np.linspace(0, 1, len(columnsMu))]
    
            colorsList = []
            #for label in labels:
            #    colorsList.append(colors[label])
            #plt.scatter(x, y, s = 80, alpha=.3, color=colorsList)
            
            
            targetVal = np.asanyarray(Sum_of_squared_distances)
            stdVal = (targetVal - targetVal.min()) / (targetVal.max() - targetVal.min())
            scaledSSE = stdVal * (max(x) - min(x)) + min(x)
            
            targetVal = np.asanyarray(K)
            stdVal = (targetVal - targetVal.min()) / (targetVal.max() - targetVal.min())
            scaledK = stdVal * (max(x) - min(x)) + min(x)
            
            axs[i,j].plot(scaledK, scaledSSE, 'kx--', lw=1, alpha=0.5)
            #plt.xlabel('k')
            #plt.ylabel('Sum_of_squared_distances')
            #plt.title('Elbow Method For Optimal k')
            #plt.show()
            
            
            
            for k in range(len(labels)):
                axs[i,j].scatter(x[k], y[k], s = 80, alpha=.4, color=colors[labels[k]])
                axs[i,j].annotate(str(k+1), (x[k], y[k]))
            
            if i == 0 and j == 0:
                axs[i,j].set_title('GD')
            if i == 0 and j == 1:
                axs[i,j].set_title('IGD')
            if i == 0 and j == 2:
                axs[i,j].set_title('HV')
    
            if i == 0 and j == 0:
                axs[i,j].set_ylabel("Morris LHS")
    
            if i == 1 and j == 0:
                axs[i,j].set_ylabel("Morris")          
    
            if i == 2 and j == 0:
                axs[i,j].set_ylabel("Sobol")          
    
            # GO NEXT PLOT            
            j = j + 1
            if j%3 == 0:
                j = 0
                i = i + 1
                
    plt.tight_layout()
    plt.savefig(os.path.join(plots, "Plot_Cluster_"+plotAlgo+"_Func.pdf"),bbox_inches='tight')
    plt.savefig(os.path.join(plots, "Plot_Cluster_"+plotAlgo+"_Func.png"),bbox_inches='tight')
    plt.show()

#%% Cluster Plot param


clusterMOEADk = [[3,5,3],
                 [4,5,4],
                 [3,5,4]]

clusterNSGAk = [[3,3,2],
                [3,3,2],
                [3,4,2]]



for algo in AlgoList:
    root = os.path.join(folder,algo)
    listdir = os.listdir(root)
    
    fig, axs = plt.subplots(3, 3, figsize=(8, 8))#, sharey=True, sharex=True)
    
    checkTest = 'var'
    i = 0
    j = 0
    for file in listdir:
        if checkTest in file:
            print(file)
            pathData = os.path.join(root,file)
            data = pd.read_csv(pathData)
            columns = data.columns.tolist()
            columnsMu = [colName for colName in columns if '.1' not in colName]
            columnsSi = [colName for colName in columns if '.1' in colName]
            dataMu = data[columnsMu]
            dataSi = data[columnsSi]
            
            dataMuT = dataMu.T       
            dataMuTnorm=(dataMuT-dataMuT.min())/(dataMuT.max()-dataMuT.min())
            dataSiT = dataSi.T       
            dataSiTnorm=(dataSiT-dataSiT.min())/(dataSiT.max()-dataSiT.min())
            
            dataMuTnorm = dataMuTnorm.fillna(0)
            dataSiTnorm = dataSiTnorm.fillna(0)
            
            dataMunorm = dataMuTnorm.T
            dataSinorm = dataSiTnorm.T
            
            dataClusterFunc  = pd.concat([dataMunorm, dataSinorm], axis=1)
            dataClusterParams = pd.concat([dataMuTnorm, dataSiTnorm], axis=0)
            
            dataCluster = dataClusterParams
            
            from sklearn.cluster import KMeans
            from sklearn.metrics import silhouette_samples, silhouette_score
            Sum_of_squared_distances = []
            silhouette_avg_values = []
            Xmat = dataCluster.to_numpy()
            K = range(2,len(dataCluster))
            for kIDX in K:
                km = KMeans(n_clusters=kIDX)
                km = km.fit(Xmat)
                #Sum_of_squared_distances.append(km.inertia_)
                cluster_labels = km.fit_predict(Xmat)
                silhouette_avg_values.append(silhouette_score(Xmat, cluster_labels))
            
            Sum_of_squared_distances = silhouette_avg_values
            np.argmax(silhouette_avg_values)
            #clusters_no = clusterK[i][j]
            clusters_no = 2+np.argmax(silhouette_avg_values)
            print(i,j,clusters_no) 

            kmeans = KMeans(n_clusters=clusters_no).fit(dataCluster)
            centroids = kmeans.cluster_centers_

            labels = kmeans.labels_
            
            from sklearn.decomposition import PCA
            pca = PCA(n_components=2)
            pca.fit(dataCluster.T)
            dataPCA =  pca.components_.T
            
            x =  dataPCA[:,0]
            y =  dataPCA[:,1]
                    
            if '_GD' in file:
                cmap = plt.cm.get_cmap('tab10') 
            elif '_IGD' in file:
                cmap = plt.cm.get_cmap('nipy_spectral') 
            else:
                cmap = plt.cm.get_cmap('Set1') 
                
            if 'MOEAD' in file:
                paramSTR = MOEAD
                plotAlgo = "MOEAD"
            else:
                paramSTR = NSGAIII
                plotAlgo = "NSGAIII"
                
            colors = [cmap(each) for each in np.linspace(0, 1, len(columnsMu))]
    
            colorsList = []
            #for label in labels:
            #    colorsList.append(colors[label])
            #plt.scatter(x, y, s = 80, alpha=.3, color=colorsList)
            
            targetVal = np.asanyarray(Sum_of_squared_distances)
            stdVal = (targetVal - targetVal.min()) / (targetVal.max() - targetVal.min())
            scaledSSE = stdVal * (max(x) - min(x)) + min(x)
            
            targetVal = np.asanyarray(K)
            stdVal = (targetVal - targetVal.min()) / (targetVal.max() - targetVal.min())
            scaledK = stdVal * (max(x) - min(x)) + min(x)
            
            #axs[i,j].plot(K, Sum_of_squared_distances, 'kx--', lw=1, alpha=0.5)
            axs[i,j].plot(scaledK, scaledSSE, 'kx--', lw=1, alpha=0.5)
            #plt.xlabel('k')
            #plt.ylabel('Sum_of_squared_distances')
            #plt.title('Elbow Method For Optimal k')
            #plt.show()            
            
            for k in range(len(labels)):
                axs[i,j].scatter(x[k], y[k], s = 80, alpha=.3, color=colors[labels[k]])
                
                param_no = str(k+1)
                if k >= len(labels)/2:
                    param_no = "i"+str(int(k-len(labels)/2)+1)
                #print(param_no)
                    
                axs[i,j].annotate(param_no, (x[k], y[k]))
                #axs[i,j].annotate(str(k+1), (x[k], y[k]))
            
            if i == 0 and j == 0:
                axs[i,j].set_title('GD')
            if i == 0 and j == 1:
                axs[i,j].set_title('IGD')
            if i == 0 and j == 2:
                axs[i,j].set_title('HV')
    
            if i == 0 and j == 0:
                axs[i,j].set_ylabel("Morris LHS")
    
            if i == 1 and j == 0:
                axs[i,j].set_ylabel("Morris")          
    
            if i == 2 and j == 0:
                axs[i,j].set_ylabel("Sobol")          
    
            # GO NEXT PLOT            
            j = j + 1
            if j%3 == 0:
                j = 0
                i = i + 1
                
    plt.tight_layout()
    plt.savefig(os.path.join(plots, "Plot_Cluster_"+plotAlgo+"_Param.pdf"),bbox_inches='tight')
    plt.savefig(os.path.join(plots, "Plot_Cluster_"+plotAlgo+"_Param.png"),bbox_inches='tight')
    plt.show()


            
#%% Statistical test of param

MOEAD = [r"$\lambda$", r"$P[\mathrm{X}]$", r"$\mathrm{X}_{\mathrm{DI}}$", r"$P[\mathrm{PM}]$", r"$\mathrm{PM}_{\mathrm{DI}}$",r"$Mode$", r"$\epsilon_N$"]
NSGAIII = [r"$\lambda$", r"$P[\mathrm{X}]$", r"$\mathrm{X}_{\mathrm{DI}}$", r"$P[\mathrm{PM}]$", r"$\mathrm{PM}_{\mathrm{DI}}$", r"$K$"]

AlgoList = ['MOEAD', 'NSGAIII']
ObjList = ['GD', 'HV', 'IGD']

for algo in AlgoList:
    root = os.path.join(folder,algo)
    listdir = os.listdir(root)
    for obj in ObjList:            
        checkTest = algo+'_'+obj+'_'
        print(checkTest)
        for file in listdir:
            if checkTest in file:
                print(file)
                pathData = os.path.join(root,file)
                data = pd.read_csv(pathData)
                columns = data.columns.tolist()
                columnsMu = [colName for colName in columns if '.1' not in colName]
                columnsSi = [colName for colName in columns if '.1' in colName]
                dataMu = data[columnsMu]
                dataSi = data[columnsSi]
                
                dataMuT = dataMu.T       
                dataMuTnorm=(dataMuT-dataMuT.min())/(dataMuT.max()-dataMuT.min())
                dataSiT = dataSi.T       
                dataSiTnorm=(dataSiT-dataSiT.min())/(dataSiT.max()-dataSiT.min())
                
                dataMuTnorm = dataMuTnorm.fillna(0)
                dataSiTnorm = dataSiTnorm.fillna(0)
                
                dataMunorm = dataMuTnorm.T
                dataSinorm = dataSiTnorm.T
                
                dataClusterFunc  = pd.concat([dataMunorm, dataSinorm], axis=1)
                dataClusterParams = pd.concat([dataMuTnorm, dataSiTnorm], axis=0)
                
                dataCluster = dataClusterParams
    
                #dataClusterParams
                dataClusterStat = dataClusterFunc
            # if check end
        #all sampling files collected            
        if algo == 'MOEAD':
            columnsStat = []
            columnsStat.extend(MOEAD) 
            columnsStat.extend(MOEAD)
        else:
            columnsStat = []
            columnsStat.extend(NSGAIII) 
            columnsStat.extend(NSGAIII)
            
        
        from scipy.stats import ttest_ind
        dataStatAll = []
        columsStatAll = []
        for col1 in dataClusterStat.columns.tolist():
            sample1  = dataClusterStat[col1].tolist()    
            dataStatCol = []
            columsStatAll = []
            i = 0
            for col2 in dataClusterStat.columns.tolist():
                sample2  = dataClusterStat[col2].tolist()
                #print(col1, col2, end='')
                tTest = ttest_ind(sample1, sample2)
                columsStatAll.append(columnsStat[i])
                columsStatAll.append(columnsStat[i])
                i = i + 1
                dataStatCol.append(tTest[0])
                dataStatCol.append(tTest[1])
            dataStatAll.append(dataStatCol)
            #print()
        
        # Create the pandas DataFrame
        dataStatAllFrame = pd.DataFrame(dataStatAll, columns = columsStatAll)
        dataStatAllFrame['index'] = columnsStat
        #dataStatAllFrame.set_index(['index'])
        saveFile = checkTest+'_Stat.csv'
        dataStatAllFrame.to_csv(saveFile)
        print('saved: ',saveFile)