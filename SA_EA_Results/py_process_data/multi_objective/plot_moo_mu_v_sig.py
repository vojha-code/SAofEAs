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

for algo in AlgoList:
    root = os.path.join(folder,algo)
    listdir = os.listdir(root)
    fig, axs = plt.subplots(3, 3, figsize=(10, 5), sharex=False, sharey=True)
    
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
            
            dataMuNormSum = dataMunorm.sum()
            dataSiNormSum = dataSinorm.sum()
            
            if '_GD' in file:
                cmap = plt.cm.get_cmap('prism') 
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
    
            if i < 2:             
                xLabel = "$\mu$"
                yLabel = "$\sigma$"
            else:
                xLabel = "$S_i$"
                yLabel = "$ST_i$"
                    
                
            colors = [cmap(each) for each in np.linspace(0, 1, len(columnsMu))]      
            colors = ['black', 'gray']
            width = 0.35
            x = dataMuNormSum.tolist()
            y = dataSiNormSum.tolist()
            x =  np.asarray(x)
            x =  (x - np.min(x)) / (np.max(x) - np.min(x))
            y =  np.asarray(y)
            y =  (y - np.min(y)) / (np.max(y) - np.min(y))
            
            xpy = x+y
            indexSort = np.argsort(xpy)
            
            paramSort = [paramSTR[idx] for idx in indexSort]
            xSort = [x[idx] for idx in indexSort]
            ySort = [y[idx] for idx in indexSort]
        
            #axs[i,j].bar(paramSTR, x, width, label =xLabel, color=colors[0], alpha=0.7)
            #axs[i,j].bar(paramSTR, y, width, bottom = x, label=yLabel, color=colors[1],alpha=0.7)
            
            axs[i,j].bar(paramSort, xSort, width, label =xLabel, color=colors[0], alpha=1)
            axs[i,j].bar(paramSort, ySort, width, bottom = xSort, label=yLabel, color=colors[1],alpha=1)
            
            if j == 2 and i == 2:
                axs[i,j].legend(ncol=1, loc='upper left')
            else:
                axs[i,j].legend(ncol=1)
    
            if i == 0 and j == 0:
                axs[i,j].set_title('GD')
            if i == 0 and j == 1:
                axs[i,j].set_title('IGD')
            if i == 0 and j == 2:
                axs[i,j].set_title('HV')
                
            if i == 0 and j == 0:
                axs[i,j].set_ylabel('Morris LHS\nScore')
            if i == 1 and j == 0:
                axs[i,j].set_ylabel('Morris\nScore')
            if i == 2 and j == 0:
                axs[i,j].set_ylabel('Sobol\nScore')
            
            # GO NEXT PLOT            
            i = i + 1
            if i%3 == 0:
                i = 0
                j = j + 1
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots, "Plot_Bar_"+plotAlgo+".pdf"),bbox_inches='tight')
    plt.savefig(os.path.join(plots, "Plot_Bar_"+plotAlgo+".png"),bbox_inches='tight')
    plt.show()

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
    plt.savefig(os.path.join(plots, "Plot_Cluster_"+plotAlgo+"_Param.pdf"),bbox_inches='tight')
    plt.savefig(os.path.join(plots, "Plot_Cluster_"+plotAlgo+"_Param.png"),bbox_inches='tight')
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
                cmap = plt.cm.get_cmap('prism') 
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
                
            plt.tight_layout()
            plt.savefig(os.path.join(plots, file+".pdf"),bbox_inches='tight')
            plt.savefig(os.path.join(plots, file+".png"),bbox_inches='tight')
            plt.show()