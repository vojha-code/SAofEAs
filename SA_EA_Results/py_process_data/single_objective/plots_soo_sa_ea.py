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
#%%
root = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_SOO\SOO_10K_FunEval_33_Funs'
listdir = os.listdir(root)

DE = [r"$\lambda$", r"$\mathbf{b}_{\mathrm{type}}$", r"$\mathbf{b}\lambda_{\mathrm{ratio}}$", r"$\mathrm{X}$", r"$P[\mathrm{X}]$", r"$\beta_{\mathrm{min}}$", r"$\beta_{\mathrm{max}}$"]
CMAES = [r"$\lambda$", r"$\mu\lambda_{\mathrm{ratio}}$", r"$\sigma_0$", r" $\alpha_{\mu}$", r"$\sigma_{0-scale}$"]

#%%
markerList = ["o", "v", "*", "X", "D", "+", "3", "s"]

#plt.style.use('seaborn-whitegrid')
fig, axs = plt.subplots(3, 2, figsize=(7, 8), sharey=True, sharex=True)

checkTest = '_Sample1K_Term10K'
i = 0
j = 0
for file in listdir:
    if checkTest in file:
        print(i,j,file)
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
        
        if 'DE' in file:
            cmap = plt.cm.get_cmap('tab10') 
        else:
            cmap = plt.cm.get_cmap('Accent') 
        colors = [cmap(each) for each in np.linspace(0, 1, len(columnsMu))]

        for paramIndex in range(len(columnsMu)):
            paramX = columnsMu[paramIndex]
            paramY = columnsSi[paramIndex]
            if 'DE' in file:
                paramSTR = DE[paramIndex]
            else:
                paramSTR = CMAES[paramIndex]
                 
            x = dataMunorm[paramX]
            y = dataSinorm[paramY]
            xmean = np.mean(dataMunorm[paramX].dropna().tolist())
            ymean = np.mean(dataSinorm[paramY].dropna().tolist())
            
            xstd = np.mean(dataMunorm[paramX].dropna().tolist())
            ystd = np.mean(dataSinorm[paramY].dropna().tolist())
            from matplotlib.patches import Ellipse
            ellipse = Ellipse(xy=(xmean, ymean), width=xstd, height=ystd, edgecolor=colors[paramIndex], fc='None', lw=1)
            axs[i,j].add_patch(ellipse)
            #axs[i,j].scatter(x, y, alpha=0.0002, s=20, marker=markerList[paramIndex], color=colors[paramIndex])
            #axs[i,j].scatter(xmean, ymean, alpha=1, s=100, marker=markerList[paramIndex], label=paramSTR, color=colors[paramIndex])

        for paramIndex in range(len(columnsMu)):
            paramX = columnsMu[paramIndex]
            paramY = columnsSi[paramIndex]
            if 'DE' in file:
                paramSTR = DE[paramIndex]
            else:
                paramSTR = CMAES[paramIndex]
                 
            x = dataMunorm[paramX]
            y = dataSinorm[paramY]
            xmean = np.mean(dataMunorm[paramX].dropna().tolist())
            ymean = np.mean(dataSinorm[paramY].dropna().tolist())
            
            #xstd = np.mean(dataMunorm[paramX].dropna().tolist())
            #ystd = np.mean(dataSinorm[paramY].dropna().tolist())
            #from matplotlib.patches import Ellipse
            #ellipse = Ellipse(xy=(xmean, ymean), width=xstd, height=ystd, edgecolor=colors[paramIndex], fc='None', lw=1)
            #axs[i,j].add_patch(ellipse)
            #axs[i,j].scatter(x, y, alpha=0.0002, s=20, marker=markerList[paramIndex], color=colors[paramIndex])
            axs[i,j].scatter(xmean, ymean, alpha=1, s=80, marker=markerList[paramIndex], label=paramSTR, color=colors[paramIndex])
            #axs[i,j].plot((0,1), (1,1), alpha=0.1, color='gray', lw=0.5)
            #axs[i,j].plot((1,1), (1,0), alpha=0.1, color='gray', lw=0.5)
            axs[i,j].set_ylim(-0.05,1.05)
            axs[i,j].set_xlim(-0.05,1.05)
            
        if i == 0 and j == 0:
            axs[i,j].set_title('CMAES')
        if i == 0 and j == 1:
            axs[i,j].set_title('DE')

        if i == 0 and j == 0:
            axs[i,j].set_ylabel("Morris LHS $\sigma$")          
        if i == 0:
            axs[i,j].set_xlabel("Morris LHS  $\mu$ ") 
        
        if i == 1 and j == 0:
            axs[i,j].set_ylabel("Morris  $\sigma$")          
        if i == 1:
            axs[i,j].set_xlabel("Morris  $\mu$") 
            
        if i == 2:
            axs[i,j].set_xlabel("Sobol  $S_i$")         
            
        if i == 2 and j == 0: 
            axs[i,j].set_ylabel("Sobol $ST_i$")          
            axs[i,j].legend(bbox_to_anchor=(1, 0.888), labelspacing=2.5, loc="upper center",  bbox_transform=fig.transFigure, ncol=1)
        
        if i == 2 and j == 1: 
            axs[i,j].legend(bbox_to_anchor=(0.99, 0.55), labelspacing=2.5, loc="upper center",  bbox_transform=fig.transFigure, ncol=1)
        
        # GO NEXT PLOT            
        j = j + 1
        if j%2 == 0:
            j = 0
            i = i + 1

plt.tight_layout()
plt.savefig(os.path.join(plots, "Plot_Scatter_CMAES-DE.pdf"),bbox_inches='tight')
plt.savefig(os.path.join(plots, "Plot_Scatter_CMAES-DE.png"),bbox_inches='tight')
plt.show()

#%% Bar plot

#plt.style.use('seaborn-whitegrid')
fig, axs = plt.subplots(2, 4, figsize=(12, 5), sharey=True)

checkTest = '_Sample1K_Term10K'

Algo = ['ES', 'DE']

i = 0
j = 0
for algo in Algo:
    xAll = []
    yAll = []
    pAll = []
    for file in listdir:
        if algo+checkTest in file:
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
            
            if 'DE' in file:
                cmap = plt.cm.get_cmap('prism') 
                paramSTR = DE
            else:
                cmap = plt.cm.get_cmap('nipy_spectral') 
                paramSTR = CMAES
    
            if j < 2:             
                xLabel = "$\mu$"
                yLabel = "$\sigma$"
            
            if j == 2 :
                xLabel = "$S_i$"
                yLabel = "$ST_i$"
                
            colors = [cmap(each) for each in np.linspace(0, 1, len(columnsMu))] 
            colors = ['black', 'lightgray']
            width = 0.35
            x = dataMuNormSum.tolist()
            y = dataSiNormSum.tolist()
            x =  np.asarray(x)
            x =  (x - np.min(x)) / (np.max(x) - np.min(x))
            y =  np.asarray(y)
            y =  (y - np.min(y)) / (np.max(y) - np.min(y))
            
            xAll.append(x)
            yAll.append(y)
            pAll.append(paramSTR)        
            
            xpy = x+y
            indexSort = np.argsort(xpy)
            
            paramSort = [paramSTR[idx] for idx in indexSort]
            xSort = [x[idx] for idx in indexSort]
            ySort = [y[idx] for idx in indexSort]
            
            #axs[i,j].bar(paramSTR, x, width, label =xLabel, color=colors[0], alpha=0.7)
            #axs[i,j].bar(paramSTR, y, width, bottom = x, label=yLabel, color=colors[1],alpha=0.7)
            
            axs[i,j].bar(paramSort, xSort, width, label =xLabel, color='w', alpha=1, edgecolor=colors[0], linewidth=1)
            axs[i,j].bar(paramSort, ySort, width, bottom = xSort, label=yLabel, color=colors[1], alpha=1, edgecolor=colors[0], linewidth=1)
    
            #if j == 0 and i == 1:
            #    axs[i,j].legend(ncol=1, loc='upper center')
            #else:
            axs[i,j].legend(ncol=1)
            
            print(i,j,algo,file)
    
            if i == 0 and j == 0:
                axs[i,j].set_title('Morris LHS')
            if i == 0 and j == 1:
                axs[i,j].set_title('Morris ')
            if i == 0 and j == 2:
                axs[i,j].set_title('Sobol ')
                
            if i == 0 and j == 0:
                axs[i,j].set_ylabel('CMAES\nScore')
            if i == 1 and j == 0:
                axs[i,j].set_ylabel('DE\nScore')
            
            # GO NEXT ALGO PLOT            
            j = j + 1
            
    # PLOT ALL    
    if j == 3 :
        xLabel = "direct"
        yLabel = "interaction"
    
    x =  np.zeros(len(xAll[0]))
    y =  np.zeros(len(yAll[0]))
    for k in range(len(xAll)):
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
    
    print(i,j,algo,file)
    
    
    axs[i,j].bar(paramSort, xSort, width, label =xLabel, color='w', alpha=1, edgecolor=colors[0], linewidth=1)
    axs[i,j].bar(paramSort, ySort, width, bottom = xSort, label=yLabel, color=colors[1], alpha=1, edgecolor=colors[0], linewidth=1)
    
    axs[i,j].legend(ncol=1)
    
    if i == 0 and j == 3:
        axs[i,j].set_title('LHS + Morris + Sobol ')        
    
    if j == 3 and i == 0:
        print('>>>>>', algo)
    # GO NEXT ALGO PLOT            
    i = i + 1
    j = 0
        


plt.tight_layout()
plt.savefig(os.path.join(plots, "Plot_Bar_CMAES-DE.pdf"),bbox_inches='tight')
plt.savefig(os.path.join(plots, "Plot_Bar_CMAES-DE.png"),bbox_inches='tight')
plt.show()
        
#%% Line plot of all data

root = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_SOO'
listdir = os.listdir(root)

checkTest = '_All.csv'
for file in listdir:
    if checkTest in file:
        print(file)
        pathData = os.path.join(root,file)
        data = pd.read_csv(pathData)    
        columns = data.columns.tolist()
        reduce = 5
        fig, axs = plt.subplots(1, len(columns)-reduce, figsize=(2*(len(columns)-reduce), (len(columns)-reduce)/2+1), sharey =True)
        #fig, axs = plt.subplots(1, len(columns)-4, sharey =True)
        
        if 'DE' in file:
            cmap = plt.cm.get_cmap('nipy_spectral') 
            paramSTR = DE
            scoreType = 'DE'
            plotName =  "DE_All"
        else:
            cmap = plt.cm.get_cmap('nipy_spectral') 
            paramSTR = CMAES
            scoreType = 'CMAES'
            plotName =  "CMA_ES_All"

        colors = [cmap(each) for each in np.linspace(0, 1, len(paramSTR))]                  
        j = 0
        for paramIndex in range(3,len(columns)-2):
            panamLabel = paramSTR[paramIndex-3]
            colorsVal = colors[paramIndex-3]
            print(paramIndex,'  ',columns[paramIndex], '=', panamLabel)
            paramX = columns[paramIndex]
            
            
            dataXY = data[[paramX, 'BestSol_mean']]
            dataXYGroup = dataXY.groupby(paramX, as_index=False).mean()
    
            x = dataXYGroup[paramX].tolist()
            #y = dataXYGroup['BestSol_mean'].tolist()
            bins = np.linspace(min(x), max(x), num=50).tolist()
            groupsMean = dataXYGroup.groupby(pd.cut(dataXYGroup[paramX], bins)).mean()
            groupsStd = dataXYGroup.groupby(pd.cut(dataXYGroup[paramX], bins)).std()
            x = groupsMean[paramX].tolist()
            #x = [k for k in range(len(groupsMean[paramX].tolist()))]
            y = groupsMean['BestSol_mean'].tolist()
            yStd = groupsStd['BestSol_mean'].tolist()
            
            #plt.scatter(x,y)
            yVal =  np.asarray(y)
            yNorm =  1 - (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))
            yVal =  np.asarray(yStd)
            yStdNorm =  1 - (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))
            axs[j].plot(x,yNorm, label=panamLabel, alpha=0.7,  color=colorsVal, linestyle="",marker="o")
            #yNorm = yNorm.tolist()
            #axs[i,j].scatter(x,yNorm)#, label=panamLabel, color=colorsVal, alpha=0.7)
            #axs[j].set_ylim(0.5,1.0)
            if 'DE' in file:
                axs[j].legend(framealpha=1, fontsize=14)
            else:
                axs[j].legend(framealpha=1)

            j = j + 1
        
        if 'DE' in file:
            axs[0].set_ylabel(scoreType+'  Score', fontsize=14)
        else:
            axs[0].set_ylabel(scoreType+'  Score')
        
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots, plotName+".pdf"),bbox_inches='tight')
        plt.savefig(os.path.join(plots, plotName+".png"),bbox_inches='tight')
        plt.show()
    
#%% Line plot of all data (Specific to methods)

root = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_SOO'
listdir = os.listdir(root)

from scipy.ndimage import gaussian_filter1d

checkTest = '_All.csv'
for file in listdir:
    if checkTest in file:
        print(file)
        pathData = os.path.join(root,file)
        data = pd.read_csv(pathData)    
        columns = data.columns.tolist()
        reduce = 5
        fig, axs = plt.subplots(1, len(columns)-reduce, figsize=(2*(len(columns)-reduce), (len(columns)-reduce)/2+1), sharey =True)
        #fig, axs = plt.subplots(1, len(columns)-4, sharey =True)
        
        if 'DE' in file:
            paramSTR = DE
            scoreType = 'DE'
            plotName =  "DE_All"
        else:
            paramSTR = CMAES
            scoreType = 'CMAES'
            plotName =  "CMA_ES_All"

        cmap = plt.cm.get_cmap('nipy_spectral') 
        #colors = [cmap(each) for each in np.linspace(0, 1, len(paramSTR))]                  
        colors = [cmap(each) for each in np.linspace(0, 1, 7)]                  
        j = 0
        for paramIndex in range(3,len(columns)-2):
            panamLabel = paramSTR[paramIndex-3]
            #colorsVal = colors[paramIndex-3]
            print(paramIndex,'  ',columns[paramIndex], '=', panamLabel)
            paramX = columns[paramIndex]
            
            
            dataXY = data[[paramX, 'BestSol_mean']]
            dataXYGroup = dataXY.groupby(paramX, as_index=False).mean()
    
            x = dataXYGroup[paramX].tolist()
            #y = dataXYGroup['BestSol_mean'].tolist()
            bins = np.linspace(min(x), max(x), num=50).tolist()
            groupsMean = dataXYGroup.groupby(pd.cut(dataXYGroup[paramX], bins)).mean()
            groupsMean = groupsMean.dropna()
            groupsStd = dataXYGroup.groupby(pd.cut(dataXYGroup[paramX], bins)).std()
            x = groupsMean[paramX].tolist()
            #x = [k for k in range(len(groupsMean[paramX].tolist()))]
            y = groupsMean['BestSol_mean'].tolist()
            yStd = groupsStd['BestSol_mean'].tolist()
            
            #plt.scatter(x,y)
            yVal =  np.asarray(y)
            yNorm =  1 - (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))
            yVal =  np.asarray(yStd)
            yStdNorm =  1 - (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))
            
            #axs[j].plot(x,yNorm, label=panamLabel, alpha=1,  color=colors[0], linestyle="",marker="o")
            
            #Paper plot
            #axs[j].plot(x,yNorm, alpha=1,  color=colors[0], linestyle="",marker="o")
            
            yNorm_smoothed = gaussian_filter1d(yNorm, sigma=2)
            axs[j].plot(x,yNorm_smoothed, alpha=1,  color=colors[0], linestyle="-",marker="o")
            
            #yNorm = yNorm.tolist()
            #axs[i,j].scatter(x,yNorm)#, label=panamLabel, color=colorsVal, alpha=0.7)
            #axs[j].set_ylim(0.5,1.0)
            #if 'DE' in file:
            #    axs[j].legend(framealpha=1, fontsize=14)
            #else:
            #    axs[j].legend(framealpha=1)
            
            #colorsVal =  [colors['k'],colors['b'],colors['r']]
            Methods = ['MorrisLHS','Morris','SOBOL']#,
            MethodLebel = ['l','m','s']
            for k in range(len(Methods)):
                dataXY = data[(data[columns[0]] == Methods[k])][[paramX, 'BestSol_mean']]
            
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
                y = groupsMean['BestSol_mean'].tolist()
                yStd = groupsStd['BestSol_mean'].tolist()
                
                #plt.scatter(x,y)
                yVal =  np.asarray(y)
                yNorm =  1 - (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal))#  reversing the value
                yVal =  np.asarray(yStd)
                yStdNorm =  1 - (yVal - np.min(yVal)) / (np.max(yVal) - np.min(yVal)) #  reversing the value
                
                #axs[i,j].plot(x,yNorm, label=panamLabel+MethodLebel[k], color=colors[k+1])
                #axs[j].plot(x,yNorm, label=panamLabel+MethodLebel[k], alpha=0.4,  color=colors[k+1], linestyle="",marker="o")
                #axs[j].legend(fontsize=16)
                
                #axs[j].plot(x,yNorm, label=panamLabel, alpha=0.4,  color=colors[k+1], linestyle="-",marker="o")
                
                #Paper plot
                #axs[j].plot(x,yNorm, alpha=0.5, color=colors[k+1])
                
                yNorm_smoothed = gaussian_filter1d(yNorm, sigma=2)
                axs[j].plot(x,yNorm_smoothed, alpha=0.5, color=colors[k+1])
                #axs[i,j].plot(x,yNorm, color=colors[k], linestyle="-",marker="o")
                #if i == 2:
                axs[j].set_xlabel(panamLabel, fontsize=16)

            j = j + 1
        
        if 'DE' in file:
            axs[0].set_ylabel(scoreType+'  Score', fontsize=14)
        else:
            axs[0].set_ylabel(scoreType+'  Score')
        
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots, plotName+".pdf"),bbox_inches='tight')
        plt.savefig(os.path.join(plots, plotName+".png"),bbox_inches='tight')
        plt.show()

#%% Cluster plot 
root = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_SOO\SOO_10K_FunEval_33_Funs'
listdir = os.listdir(root)

markerList = ["o", "v", "*", "X", "D", "+", "3", "s"]



clusterK = [[4,4],
            [5,4],
            [6,5]]



#plt.style.use('seaborn-whitegrid')
fig, axs = plt.subplots(3, 2, figsize=(7, 8))#, sharey=True, sharex=True)

checkTest = '_Sample1K_Term10K'
i = 0
j = 0
# for index in range(6):
#     print(i,j,clusterK[i][j])     
#     # GO NEXT PLOT            
#     j = j + 1
#     if j%2 == 0:
#         j = 0
#         i = i + 1
            
for file in listdir:
    if checkTest in file:
        #print(file)
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
                
        if 'DE' in file:
            cmap = plt.cm.get_cmap('tab10') 
        else:
            cmap = plt.cm.get_cmap('Accent') 
        colors = [cmap(each) for each in np.linspace(0, 1, 10)]

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
        
        for k in range(len(labels)):
            axs[i,j].scatter(x[k], y[k], s = 80, alpha=.5, color=colors[labels[k]])
            
            param_no = str(k+1)
            #if k > len(labels)/2:
            #    param_no = str(int(k-len(labels)/2)+1)
                            
            #if k < len(labels)/2:
            axs[i,j].annotate(param_no, (x[k], y[k]))
        
        if i == 0 and j == 0:
            axs[i,j].set_title('CMAES')
        if i == 0 and j == 1:
            axs[i,j].set_title('DE')

        if i == 0 and j == 0:
            axs[i,j].set_ylabel("Morris LHS")

        if i == 1 and j == 0:
            axs[i,j].set_ylabel("Morris")          

        if i == 2 and j == 0:
            axs[i,j].set_ylabel("Sobol")          

        # GO NEXT PLOT            
        j = j + 1
        if j%2 == 0:
            j = 0
            i = i + 1

plt.tight_layout()
plt.savefig(os.path.join(plots, "Plot_Cluster_CMAES-DE_Func.pdf"),bbox_inches='tight')
plt.savefig(os.path.join(plots, "Plot_Cluster_CMAES-DE_Func.png"),bbox_inches='tight')
plt.show()


#%% Cluster Plot param
root = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_SOO\SOO_10K_FunEval_33_Funs'
listdir = os.listdir(root)
#listdir.remove(listdir[1])
markerList = ["o", "v", "*", "X", "D", "+", "3", "s"]


clusterK = [[4,4],
            [5,3],
            [5,2]]


#plt.style.use('seaborn-whitegrid')
fig, axs = plt.subplots(3, 2, figsize=(7, 8))#, sharey=True, sharex=True)

checkTest = '_Sample1K_Term10K'
i = 0
j = 0
for file in listdir:
    if checkTest in file:
        print(i,j,file)
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
                
        if 'DE' in file:
            cmap = plt.cm.get_cmap('tab10') 
        else:
            cmap = plt.cm.get_cmap('Accent') 
        colors = [cmap(each) for each in np.linspace(0, 1, 10)]

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
        
        for k in range(len(labels)):
            axs[i,j].scatter(x[k], y[k], s = 80, alpha=.5, color=colors[labels[k]])
            
            param_no = str(k+1)
            if k >= len(labels)/2:
                param_no = "i"+str(int(k-len(labels)/2)+1)
                
            axs[i,j].annotate(param_no, (x[k], y[k]))
        
        if i == 0 and j == 0:
            axs[i,j].set_title('CMAES')
        if i == 0 and j == 1:
            axs[i,j].set_title('DE')

        if i == 0 and j == 0:
            axs[i,j].set_ylabel("Morris LHS")

        if i == 1 and j == 0:
            axs[i,j].set_ylabel("Morris")

        if i == 2 and j == 0:
            axs[i,j].set_ylabel("Sobol")

        # GO NEXT PLOT            
        j = j + 1
        if j%2 == 0:
            j = 0
            i = i + 1

plt.tight_layout()
plt.savefig(os.path.join(plots, "Plot_Cluster_CMAES-DE_Param.pdf"),bbox_inches='tight')
plt.savefig(os.path.join(plots, "Plot_Cluster_CMAES-DE_Param.png"),bbox_inches='tight')
plt.show()

#%% Statistical test of param
root = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_SOO\SOO_10K_FunEval_33_Funs'
listdir = os.listdir(root)

DE = [r"$\lambda$", r"$\mathbf{b}_{\mathrm{type}}$", r"$\mathbf{b}\lambda_{\mathrm{ratio}}$", r"$\mathrm{X}$", r"$P[\mathrm{X}]$", r"$\beta_{\mathrm{min}}$", r"$\beta_{\mathrm{max}}$"]
CMAES = [r"$\lambda$", r"$\mu\lambda_{\mathrm{ratio}}$", r"$\sigma_0$", r" $\alpha_{\mu}$", r"$\sigma_{0-scale}$"]


for algo in ['ES','DE']:
    checkTest = algo+'_Sample1K_Term10K'
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
    if algo == 'DE':
        saveFile = 'DE_Stat.csv'
        columnsStat = []
        columnsStat.extend(DE)
        columnsStat.extend(DE)
    else:
        saveFile = 'CMAES_Stat.csv'
        columnsStat = []
        columnsStat.extend(CMAES)
        columnsStat.extend(CMAES)
        
    
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
    dataStatAllFrame.to_csv(saveFile)
    print('saved: ',saveFile)
        
        