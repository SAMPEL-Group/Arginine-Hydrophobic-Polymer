import pandas as pd
import numpy as np
import seaborn as sns
import hdbscan
from numpy import savetxt
import MDAnalysis as md
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.font_manager as font_manager
import matplotlib.ticker as ticker
import copy

font = "Helvetica"
hfont = {'fontname':'Helvetica'}
big = 36
medium = 28
small = 12
xsmall = 10
ms=12
lw=2

normalizeData = False

nW=12
nR=3

titles=["93lys", "185lys", "93glu", "185glu", "94arg-lys", "94arg-glu", "94lys-glu", "186arg-lys", "186arg-glu", "186lys-glu"]

#%% Import Data

for title in titles:
    windows=np.arange(0,nW,1)
    runs=np.arange(1,nR+1,1)
    
    def openFile(title, window, run):
        data = pd.read_table("data/polystat-{}-reus-{}-{}.xvg".format(title, window, run), skipinitialspace=True, index_col=False, skiprows=28, sep=" ", names=['time','ee','rg','eig1','eig2','eig3'])
        return data
    
    catData=[]
    catID=[]
    for run in runs:
        for window in windows:
            data = openFile(title, window, run)
            keepTrack = np.zeros_like(data.index)
            keepTrack = np.expand_dims(keepTrack, axis=1)
            keepTrack = np.concatenate((keepTrack, keepTrack), axis=1)
            keepTrack[:,0] = window
            keepTrack[:,1] = run
            keepTrack = pd.DataFrame(keepTrack, columns=["window", "run"])
            
            catData.append(data)
            catID.append(keepTrack)
    
    tracker = pd.concat(catID, axis=0)
    data = pd.concat(catData, axis=0)
    origData = copy.deepcopy(data)
    
    if normalizeData == True:
        data.rg = (data.rg-np.min(data.rg))/np.max((data.rg-np.min(data.rg)))
        data.ee = (data.ee-np.min(data.ee))/np.max((data.ee-np.min(data.ee)))
    
    inputData = np.vstack((data.eig1, data.eig2, data.eig3)).T
    minSize = int(0.01*len(inputData[:,0]))
    
    #%% Cluster Data
    clusterer = hdbscan.HDBSCAN(algorithm='best', alpha=1.0, approx_min_span_tree=True, 
                                gen_min_span_tree=True, metric='euclidean', 
                                min_cluster_size=100, min_samples=50, p=None,
                                cluster_selection_method='leaf', allow_single_cluster=True)
    
    clusterer.fit(inputData)
    palette = sns.color_palette(n_colors=100)
    
    cluster_colors = [sns.desaturate(palette[col], sat)
                      if col >= 0 else (0.5, 0.5, 0.5) for col, sat in
                      zip(clusterer.labels_, clusterer.probabilities_)]
    
    cluster_colors = [x for x in cluster_colors if x != (0.5, 0.5, 0.5)] # Drop noise
    
    labels = clusterer.labels_
    
    hairpins = origData[labels == 2].index
    plotData = data[labels >= 0]
    
    fig, ax = plt.subplots(figsize=(9,6))
    
    ax.scatter(plotData.rg, plotData.ee, c=cluster_colors)
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    plt.xticks(fontsize=medium, **hfont)
    plt.yticks(fontsize=medium, **hfont)
    plt.xlabel('$R_g$ [nm]', **hfont, fontsize=big)
    plt.ylabel("End-to-End [nm]", **hfont, fontsize=big)
    plt.ylim(0,3.25)
    plt.xlim(0.35, 0.925)
    ax.tick_params(axis='both', which='both', direction='in', labelsize=medium)
    ax.tick_params(axis='both', which='major', length=6, width=2)
    ax.tick_params(axis='both', which='minor', length=4, width=2)
    plt.tight_layout()
    plt.savefig('figures/HDBSCAN-{}.png'.format(title), dpi=300)
    plt.show()
    
    labels = clusterer.labels_
    probability = clusterer.probabilities_
    clusterProbabilities = pd.DataFrame(np.vstack((data.index, data.time, labels, probability)).T, columns=["frame", "time", "cluster", "prob"])
    highProbs = clusterProbabilities[clusterProbabilities.prob == 1].sort_values('cluster')
    uniqueClusters = np.unique(highProbs.cluster)
    
    frames=[]
    times=[]
    for x in uniqueClusters:
        frames.append(int(np.asarray(highProbs[highProbs.cluster == x].index)[0]))
        times.append(int(np.asarray(highProbs[highProbs.cluster == x].time)[0]))
    
    print("Frames for cluster centers are at: {}".format(frames))
    
    from numpy import savetxt
    header = "frame cluster probability"
    np.savetxt("output/{}-hdbscan.dat".format(title), clusterProbabilities, header=header, fmt='%.3f', delimiter=' ', comments='')
    
    extractor = np.asarray(tracker)
    
    frameOutput=[]
    for i in range(len(frames)):
        x = frames[i]
        time = times[i]
        window, run = extractor[x]
        print(window, run, time)
        frameOutput.append([time, window, run, i, 0])
        
    # header = "cluster window run frame"
    np.savetxt("output/{}-clusterCenters.dat".format(title), frameOutput, fmt='%.0f', delimiter=' ', comments='')
    
    # u = md.Universe(gro, traj)
    # i=0
    # for frame in frames:
    #     u.trajectory[frame]
    #     outputFrame = u.select_atoms("all")
    #     outputFrame.write("output/{}-{}-{}".format(title, int(uniqueClusters[i]), frame))
    #     i+=1
        
    # with md.Writer("output/{}-clusters.pdb".format(title), multiframe=True) as W:
    #     for frame in frames:
    #         u.trajectory[frame] # point to desired frame
    #         outputFrame = u.select_atoms("all")
    #         W.write(outputFrame)
