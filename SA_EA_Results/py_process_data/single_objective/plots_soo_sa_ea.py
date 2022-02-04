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
root = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_SOO\de-cma-es-10K_nfeval'
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
fig, axs = plt.subplots(2, 3, figsize=(10, 4), sharey=True)

checkTest = '_Sample1K_Term10K'
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
        
        if 'DE' in file:
            cmap = plt.cm.get_cmap('prism') 
            paramSTR = DE
        else:
            cmap = plt.cm.get_cmap('nipy_spectral') 
            paramSTR = CMAES

        if j < 2:             
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
        
        #if j == 0 and i == 1:
        #    axs[i,j].legend(ncol=1, loc='upper center')
        #else:
        axs[i,j].legend(ncol=1)

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
        
        # GO NEXT PLOT            
        i = i + 1
        if i%2 == 0:
            i = 0
            j = j + 1


plt.tight_layout()
plt.savefig(os.path.join(plots, "Plot_Bar_CMAES-DE.pdf"),bbox_inches='tight')
plt.savefig(os.path.join(plots, "Plot_Bar_CMAES-DE.png"),bbox_inches='tight')
plt.show()
        
#%% Line plot of all data

root = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_SOO\de-cma-es-10K_nfeval'
listdir = os.listdir(root)

checkTest = '_All.csv'
for file in listdir:
    if checkTest in file:
        print(file)
        pathData = os.path.join(root,file)
        data = pd.read_csv(pathData)    
        columns = data.columns.tolist()
        fig, axs = plt.subplots(1, len(columns)-4, figsize=(2*(len(columns)-4), (len(columns)-4)/2+1), sharey =True)
        #fig, axs = plt.subplots(1, len(columns)-4, sharey =True)
        
        if 'DE' in file:
            cmap = plt.cm.get_cmap('prism') 
            paramSTR = DE
            scoreType = 'DE'
        else:
            cmap = plt.cm.get_cmap('nipy_spectral') 
            paramSTR = CMAES
            scoreType = 'CMAES'

        colors = [cmap(each) for each in np.linspace(0, 1, len(paramSTR))]                  
        j = 0
        for paramIndex in range(3,len(columns)-1):
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
            axs[j].plot(x,yNorm, label=panamLabel, color=colorsVal, linestyle="",marker="o")
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
        plt.savefig(os.path.join(plots, file+".pdf"),bbox_inches='tight')
        plt.savefig(os.path.join(plots, file+".png"),bbox_inches='tight')
        plt.show()
    
#%% Cluster plot 
root = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_SOO\de-cma-es-10K_nfeval'
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
            axs[i,j].annotate(str(k+1), (x[k], y[k]))
        
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
root = r'C:\Users\yl918888\Desktop\Evolutionary_Algorithms_Sensitivity_Analysis\SA_EA_Results\Results_SOO\de-cma-es-10K_nfeval'
listdir = os.listdir(root)

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
            axs[i,j].annotate(str(k+1), (x[k], y[k]))
        
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

        