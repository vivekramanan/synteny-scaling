"""
Post Analysis 2
Analysis on HC / KNN over a cluster variable 
Creates Figure 3C/D and Figure 4D/E
"""

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
from sklearn.metrics import jaccard_score
import matplotlib.pyplot as plt
from skbio.tree import TreeNode
from skbio import io
import networkx as nx
import plotly.graph_objects as go
from sklearn.preprocessing import StandardScaler,MinMaxScaler
import scipy.stats as stats
import seaborn as sns
from sklearn.neighbors import kneighbors_graph
import dendrogram_from_girvan_newman as dgn
from sklearn.cluster import AffinityPropagation
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.spatial import procrustes
from skbio.stats.ordination import pcoa
import plotly.graph_objects as go
import plotly.express as px
import dendrogram_from_girvan_newman as dgn
from sklearn.metrics.cluster import contingency_matrix
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, rand_score
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import squareform
from sklearn.metrics import calinski_harabasz_score
import markov_clustering as mc


def fixSixteen(df):
    """
    Fixes the 16S matrix
    @param df: 16S filename
    @return z: numpy matrix
    """
    z = df.to_numpy()[:,1:]
    z = z/100
    for i in range(len(z)):
        z[i][i] = 1.0
    
    return z 

def fixANI(ani, spec):
    """
    Fixes the ANI matrix
    @param ani: ani pandas DF
    @param spec: list of species total from 16S
    @return a: numpy matrix
    """
    ani.set_index("Unnamed: 0", inplace=True)

    ani_new = np.zeros([len(spec),len(spec)])
    for i in range(len(ani.index)):
        temp = ani.iloc[i]
        for v in range(len(temp.values)):
            if temp.values[v] != 0.0:
                a1 = spec.index(ani.columns[v])
                a2 = spec.index(ani.columns[i])
                ani_new[a1][a2] = temp.values[v]
                ani_new[a2][a1] = temp.values[v]

    ani_col = pd.DataFrame(ani_new, index=spec, columns = spec)
    a = ani_col.to_numpy()
    a = a/100

    return a

def fixCore(coreDF, spec): 
    """
    Fixes the core genes synteny metric into numpy
    @param coreDF: core pandas DF
    @param spec: list of psecies total from 16S
    @return c: numpy matrix
    """
    c = np.zeros((len(spec), len(spec)))
    for j in range(len(spec)):
        c[j][j] = 1.0
    for i in range(len(coreDF.index)):
        temp = coreDF.iloc[i]
        if temp['Species1'] in spec and temp['Species2'] in spec:
            index1 = spec.index(temp['Species1'])
            index2 = spec.index(temp['Species2'])
            if float(temp['Synteny Similarity']) >= 0.0: 
                c[index1][index2] = float(temp['Synteny Similarity'])
                c[index2][index1] = float(temp['Synteny Similarity'])
    
    return c 

def scaling(z, a, c): 
    """
    Returns the synteny scaled metric (z_scaled), the removed version (z_removed),
        and the thresholded version (z_thresh
        )
    """
    z_scaled = MinMaxScaler().fit_transform(np.dot(z, np.cov(c)))

    z_removed = np.zeros([len(spec), len(spec)])
    z_thresh = np.zeros([len(spec), len(spec)])

    for i in range(len(spec)):
        for j in range(len(spec)):
            if a[i][j] != 0.0:
                z_removed[i][j] = z_scaled[i][j]
                
            if z_scaled[i][j] > 0.82:
                z_thresh[i][j] = z_scaled[i][j]

    print("ANI Sparsity: ", 1.0 - ( np.count_nonzero(a) / float(a.size) ))
    print("Removed Sparsity: ", 1.0 - ( np.count_nonzero(z_removed) / float(a.size) ))
    print("Threshold 82 Sparsity: ", 1.0 - ( np.count_nonzero(z_thresh) / float(a.size) ))
    print("Synteny Sparsity: ", 1.0 - ( np.count_nonzero(c) / float(c.size) ))
    print("Synteny 16S Sparsity: ", 1.0 - ( np.count_nonzero(z_scaled) / float(z_scaled.size) ))
    print("16S Sparsity: ", 1.0 - ( np.count_nonzero(z) / float(z.size) ))

    return z_scaled, z_removed, z_thresh

def readInFiles(): 
    """
    Reads in the 3 main files for analysis: 16S, synteny core metric, and ANI metric
    @param z: 16S numpy matrix
    @param a: ani numpy matrix
    @param c: core synteny numpy matrix
    """
    sixteenSFile = "" # INSERT FILENAME
    coreFile = "" # INSERT FILENAME
    aniFile = "" # INSERT FILENAME

    df = pd.read_csv(sixteenSFile)
    coreDF = pd.read_csv(coreFile)
    aniDF = pd.read_csv(aniFile)

    spec = df['Unnamed: 0'].values.tolist()
    z = fixSixteen(df)
    a = fixANI(aniDF, spec)
    c = fixCore(coreDF, spec)

    return z,a,c

def hc(X,n):
    """
    @param X: distance matrix
    @param n: cluster cutoff integer
    @return clusters: cluster groupings
    @return Z: linkage matrix
    """

    Z = linkage(X, 'complete')
    clusters = cut_tree(Z, n_clusters=n)
    
    return clusters, Z

def hcOverall(z, a, c, spec, z_scaled, z_thresh, z_removed): 
    """
    Runs hierarchical clustering functions
    @param in order: 16s, ani, core synteny, list of species, 
        scaled synteny, thresholded synteny, removed synteny
    Produces two plots: hierarchical clustering silhouette scores 
        and rand scores 
    """

    silhouette_16s = []
    silhouette_syn = []
    silhouette_ani = []
    silhouette_rem = []
    silhouette_thr = []
    silhouette_onl = []

    cophenet_16s = []
    cophenet_syn = []
    cophenet_ani = []
    cophenet_rem = []
    cophenet_thr = []
    cophenet_onl = []

    ar_1 = []
    ar_2 = []
    ar_3 = []
    ar_4 = []

    maxL = 100

    for i in range(2,maxL):

        # 16s alone
        clust_16, Z_16 = hc(z,i)
        c1 = [j[0] for j in clust_16]

        # 16s + Synteny 
        clust_16syn, Z_syn = hc(z_scaled, i)
        c2 = [j[0] for j in clust_16syn]
        
        # ani 
        clust_ani, Z_ani = hc(a, i)
        c3 = [j[0] for j in clust_ani]
        
        # synteny reduced
        clust_rem, Z_rem = hc(z_removed, i)
        c4 = [j[0] for j in clust_rem]
        
        # synteny thresholded
        clust_thr, Z_thr = hc(z_thresh, i)
        c5 = [j[0] for j in clust_thr]
        
        # synteny ALONE
        clust_onl, Z_onl = hc(c, i)
        c6 = [j[0] for j in clust_onl]

        silhouette_16s.append(silhouette_score(z, c1))
        silhouette_syn.append(silhouette_score(z_scaled, c2))
        silhouette_ani.append(silhouette_score(a, c3))
        silhouette_rem.append(silhouette_score(z_removed, c4))
        silhouette_thr.append(silhouette_score(z_thresh, c5))
        silhouette_onl.append(silhouette_score(c, c6))
        
        ar_1.append(rand_score(c1, c3))
        ar_2.append(rand_score(c2, c3))
        ar_3.append(rand_score(c4, c3))
        ar_4.append(rand_score(c5, c3))

    hcPlot(silhouette_16s, silhouette_syn, silhouette_thr, silhouette_rem, silhouette_ani)
    randPlot(ar_1, ar_2, ar_3, ar_4)

def hcPlot(silhouette_16s, silhouette_syn, silhouette_thr, silhouette_rem, silhouette_ani): 
    """
    Produces the silhouette score plot
    """
    f,axes = plt.subplots(1,1, figsize=(7,4))
    plt.rcParams.update({'font.size': 30})

    x = [i for i in range(2,maxL)]

    sns.set_theme('paper')
    sns.set_style("white")

    sns.lineplot(x=x, y=silhouette_thr, label="Thr", linewidth=2, color="#FF595E")
    axes.fill_between(x, silhouette_rem, silhouette_thr, alpha=.3, color="#FF595E")

    sns.lineplot(x=x, y=silhouette_rem, label="Rem", linewidth=2, color="#3E7B4A")
    axes.fill_between(x, silhouette_syn, silhouette_rem, alpha=.3, color="#3E7B4A")

    sns.lineplot(x=x, y=silhouette_syn, label="Syn+16S", linewidth=2, color="#8DA1B9")
    axes.fill_between(x, silhouette_16s, silhouette_syn, alpha=.3, color="#8DA1B9")

    sns.lineplot(x=x, y=silhouette_16s, label="16S", linewidth=2, color="#156064")
    axes.fill_between(x, silhouette_ani, silhouette_16s, alpha=.3, color="#156064")

    sns.lineplot(x=x, y=silhouette_ani, label="ANI", linewidth=2, color="#331832")
    axes.fill_between(x, 0.15, silhouette_ani, alpha=.3, color="#331832")

    sns.despine() 
    sns.move_legend(axes, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig("silhouetteScores_complete.png", dpi=700, bbox_inches="tight") 

    print("Silhouette Scores: Descriptive Stats")
    print("Mean, STD, Median")
    print("%0.3f, %0.3f, %0.3f" %(stats(silhouette_thr)))
    print("%0.3f, %0.3f, %0.3f" %(stats(silhouette_rem)))
    print("%0.3f, %0.3f, %0.3f" %(stats(silhouette_syn)))
    print("%0.3f, %0.3f, %0.3f" %(stats(silhouette_16s)))
    print("%0.3f, %0.3f, %0.3f" %(stats(silhouette_ani)))
    print()

def randPlot(ar_1, ar_2, ar_3, ar_4):
    """
    Produces rand score box plot
    """
    x = [i for i in range(2,maxL)]

    sns.set_theme('paper')
    sns.set_style("white")
    f,axes = plt.subplots(1,1, figsize=(4,4))

    boxData = ar_1 + ar_2 + ar_3 + ar_4
    boxCols = ["16S-ANI" for i in range(98)] + ["Syn-ANI" for i in range(98)] + ["Rem-ANI" for i in range(98)] + ["Thr-ANI" for i in range(98)]
    data = pd.DataFrame()
    data['Rand Score'] = boxData
    data['Comparison'] = boxCols

    sns.stripplot(data=data, x="Comparison", y="Rand Score", hue=x+x+x+x)
    sns.despine()
    sns.move_legend(axes, "upper left", bbox_to_anchor=(1.5, 1))
    plt.savefig("randScores.png", dpi=700, bbox_inches="tight")

    print("Rand Scores: Descriptive Stats")
    print("Mean, STD, Median")
    print("%0.3f, %0.3f, %0.3f" %(stats(ar_1)))
    print("%0.3f, %0.3f, %0.3f" %(stats(ar_2)))
    print("%0.3f, %0.3f, %0.3f" %(stats(ar_3)))
    print("%0.3f, %0.3f, %0.3f" %(stats(ar_4)))
    print()

def stats(l): 
    # Descriptive statistics
    return np.mean(l), np.std(l), np.median(l)   

def weightedJac(G1, G):
    # Weighted Jaccard
    G1_dict = edgeIntoDict(G1.edges)
    G_dict = edgeIntoDict(G.edges)
    intersection = list(set(G1.edges()).intersection(set(G.edges())))
    num = 0.0
    denom = 0.0
    for l in intersection:
        num += min(G1_dict[l], G_dict[l])
        denom += max(G1_dict[l], G_dict[l])
    return num/denom

def edgeIntoDict(G_edges):
    # Turns list of edges into a dictionary
    d = {}
    edges = list(G_edges(data=True))
    for i in edges:
        d[(i[0], i[1])] = i[2]['weight']
    return d

def KNNGraphs(z, a, c, spec, z_scaled, z_thresh, z_removed):
    """
    Runs the KNN graphs section (same parameters as hcOverall)
    Produces two plots: modularity coefficient and weighted jaccard boxplot
    """
    
    modularity_16s = []
    modularity_syn = []
    modularity_ani = []
    modularity_rem = []
    modularity_thr = []

    maxL = 25

    wj_1 = []
    wj_2 = []
    wj_3 = []
    wj_4 = []

    for i in range(1,maxL):

        # 16s alone
        A1 = kneighbors_graph(z, i, mode='distance')
        G1 = nx.from_numpy_array(MinMaxScaler().fit_transform(A1.toarray()))
        c1 = nx.community.greedy_modularity_communities(G1)

        # 16s + Synteny 
        A2 = kneighbors_graph(z_scaled, i, mode='distance')
        G2 = nx.from_numpy_array(MinMaxScaler().fit_transform(A2.toarray()))
        c2 = nx.community.greedy_modularity_communities(G2)
        
        # ani 
        A3 = kneighbors_graph(a, i, mode='distance')
        G3 = nx.from_numpy_array(MinMaxScaler().fit_transform(A3.toarray()))
        c3 = nx.community.greedy_modularity_communities(G3)
        
        # synteny removed
        A4 = kneighbors_graph(z_removed, i, mode='distance')
        G4 = nx.from_numpy_array(MinMaxScaler().fit_transform(A4.toarray()))
        c4 = nx.community.greedy_modularity_communities(G4)
        
        # synteny threshold
        A5 = kneighbors_graph(z_thresh, i, mode='distance')
        G5 = nx.from_numpy_array(MinMaxScaler().fit_transform(A5.toarray()))
        c5 = nx.community.greedy_modularity_communities(G5)

        modularity_16s.append(nx.community.modularity(G1, c1))
        modularity_syn.append(nx.community.modularity(G2, c2))
        modularity_ani.append(nx.community.modularity(G3, c3))
        modularity_rem.append(nx.community.modularity(G4, c4))
        modularity_thr.append(nx.community.modularity(G5, c5))
        
        wj_1.append(weightedJac(G1, G3))
        wj_2.append(weightedJac(G2, G3))
        wj_3.append(weightedJac(G4, G3))
        wj_4.append(weightedJac(G5, G3))
    
    knnPlot(modularity_syn, modularity_thr, modularity_16s, modularity_ani, modularity_rem)
    wjPlot(wj_1, wj_2, wj_3, wj_4)


def knnPlot(modularity_syn, modularity_thr, modularity_16s, modularity_ani, modularity_rem):
    """
    Produces modularity plot
    """
    f,axes = plt.subplots(1,1, figsize=(6,5))
    plt.rcParams.update({'font.size': 30})

    x = [i for i in range(1,maxL)]

    sns.set_theme('paper')
    sns.set_style("white")

    sns.lineplot(x=x, y=modularity_thr, label="Thr", linewidth=2, color="#FF595E")
    axes.fill_between(x, modularity_syn, modularity_thr, alpha=.3, color="#FF595E")

    sns.lineplot(x=x, y=modularity_syn, label="Syn", linewidth=2, color="#8DA1B9")
    axes.fill_between(x, modularity_syn, modularity_16s, alpha=.3, color="#8DA1B9")

    sns.lineplot(x=x, y=modularity_16s, label="16S", linewidth=2, color="#156064")
    axes.fill_between(x, modularity_ani, modularity_16s, alpha=.3, color="#156064")

    sns.lineplot(x=x, y=modularity_ani, label="ANI", linewidth=2, color="#331832")
    axes.fill_between(x, modularity_rem, modularity_ani, alpha=.3, color="#331832")

    sns.lineplot(x=x, y=modularity_rem, label="Rem", linewidth=2, color="#3E7B4A")
    axes.fill_between(x, 0.0, modularity_rem, alpha=.3, color="#3E7B4A")

    sns.move_legend(axes, "upper left", bbox_to_anchor=(1.5, 1))
    sns.despine()

    plt.savefig("modularity_greedy.png", dpi=700, bbox_inches="tight")

    print("Modularity Stats: ")
    print("Mean, STD, Median")
    print("%0.3f, %0.3f, %0.3f" %(stats(modularity_thr)))
    print("%0.3f, %0.3f, %0.3f" %(stats(modularity_rem)))
    print("%0.3f, %0.3f, %0.3f" %(stats(modularity_syn)))
    print("%0.3f, %0.3f, %0.3f" %(stats(modularity_16s)))
    print("%0.3f, %0.3f, %0.3f" %(stats(modularity_ani)))
    print()

def wjPlot(wj_1, wj_2, wj_3, wj_4):
    """
    Produces weighted jaccard plot
    """
    f,axes = plt.subplots(1,1, figsize=(4,4))
    plt.rcParams.update({'font.size': 30})

    x = [i for i in range(1,maxL)]

    sns.set_theme('paper')
    sns.set_style("white")

    boxData = wj_1 + wj_2 + wj_3 + wj_4
    boxCols = ["16S-ANI" for i in range(len(x))] + ["Syn-ANI" for i in range(len(x))] + ["Rem-ANI" for i in range(len(x))] + ["Thr-ANI" for i in range(len(x))]
    data = pd.DataFrame()
    data['Weighted Jaccard'] = boxData
    data['Comparison'] = boxCols


    sns.stripplot(data=data, x="Comparison", y="Weighted Jaccard", hue=x+x+x+x)
    sns.despine()
    sns.move_legend(axes, "upper left", bbox_to_anchor=(1.5, 1))

    plt.savefig("weightedJaccard.png", dpi=700, bbox_inches="tight")

    print("Weighted Jaccard Stats: ")
    print("Mean, STD, Median")
    print("%0.3f, %0.3f, %0.3f" %(stats(wj_1)))
    print("%0.3f, %0.3f, %0.3f" %(stats(wj_2)))
    print("%0.3f, %0.3f, %0.3f" %(stats(wj_3)))
    print("%0.3f, %0.3f, %0.3f" %(stats(wj_4)))
    print()

def funcAdj(spec, df):
    """
    Creates function adjacency matrices
    @param spec: list of species
    @param df: df of choice
    @return the adjacency numpy matrix
    """
    z = np.zeros((len(spec), len(spec)))

    for i in range(len(df.index)):
        temp = df.iloc[i]
        if temp['Species1_Name'] in spec and temp['Species2_Name'] in spec:
            index1 = spec.index(temp['Species1_Name'])
            index2 = spec.index(temp['Species2_Name'])
            z[index1][index2] = temp['CosSim']
            z[index2][index1] = temp['CosSim']
    return z

def functionalBoxPlots(mgedf, metdf, argdf, virdf, coreDF):
    """
    @param dataframes of the various functions
    """
    boxData = argdf['CosSim'].values.tolist() + virdf['CosSim'].values.tolist() + mgedf['CosSim'].values.tolist() + metdf['CosSim'].values.tolist() + coreDF['Synteny Similarity'].values.tolist()
    boxCols = ['ARG' for i in range(len(argdf['CosSim']))] + ['VIR' for i in range(len(virdf['CosSim']))] + ['MGE' for i in range(len(mgedf['CosSim']))] + ['MET' for i in range(len(metdf['CosSim']))] + ['SYN' for i in range(len(coreDF['Synteny Similarity']))]

    df = pd.DataFrame()
    df['Synteny Value'] = boxData
    df['Synteny Value'] = pd.to_numeric(boxData)
    df['Comparison'] = boxCols

    sns.set_theme('paper')
    sns.set_style("white")

    f,axes = plt.subplots(1,1, figsize=(5,7))
    plt.rcParams.update({'font.size': 30})

    sns.boxplot(data=df, x='Comparison', y='Synteny Value')
    sns.despine()

    plt.savefig("synBoxPlot.png", bbox_inches="tight", dpi=700)

def functions(): 
    """
    Runs functional dataframes
    """

    mgeFile = "" # INSERT FILENAME
    metFile = "" # INSERT FILENAME
    virFile = "" # INSERT FILENAME
    argFile = "" # INSERT FILENAME

    mgedf = pd.read_csv(mgeFile, names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])
    mgedf = mgedf[1:]
    mge = funcAdj(spec, mgedf)

    metdf = pd.read_csv(metFile, names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])
    metdf = metdf[1:]
    met = funcAdj(spec, metdf)

    virdf = pd.read_csv(virFile, names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])
    virdf = virdf[1:]
    vir = funcAdj(spec, virdf)

    argdf = pd.read_csv(argFile, names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])
    argdf = argdf[1:]
    arg = funcAdj(spec, argdf)

    functionalBoxPlots(mgedf, metdf, argdf, virdf, coreDF)

    arg_scaled = MinMaxScaler().fit_transform(np.dot(z, np.cov(arg)))
    mge_scaled = MinMaxScaler().fit_transform(np.dot(z, np.cov(mge)))
    met_scaled = MinMaxScaler().fit_transform(np.dot(z, np.cov(met)))
    vir_scaled = MinMaxScaler().fit_transform(np.dot(z, np.cov(vir)))

    arg_thresh = np.zeros([len(spec), len(spec)])
    mge_thresh = np.zeros([len(spec), len(spec)])
    met_thresh = np.zeros([len(spec), len(spec)])
    vir_thresh = np.zeros([len(spec), len(spec)])

    for i in range(len(spec)):
        for j in range(len(spec)):            
            if arg_scaled[i][j] > 0.82:
                arg_thresh[i][j] = arg_scaled[i][j]
            if mge_scaled[i][j] > 0.82:
                mge_thresh[i][j] = mge_scaled[i][j]
            if met_scaled[i][j] > 0.82:
                met_thresh[i][j] = met_scaled[i][j]
            if vir_scaled[i][j] > 0.82:
                vir_thresh[i][j] = vir_scaled[i][j]

    print("Sparsities of Functional DataFrames: ")
    print("Arg Sparsity: ", 1.0 - ( np.count_nonzero(arg_thresh) / float(a.size) ))
    print("Mge Sparsity: ", 1.0 - ( np.count_nonzero(mge_thresh) / float(a.size) ))
    print("Met Sparsity: ", 1.0 - ( np.count_nonzero(met_thresh) / float(a.size) ))
    print("Vir Sparsity: ", 1.0 - ( np.count_nonzero(vir_thresh) / float(c.size) ))
    print()

    hcFunctions(arg_thresh, mge_thresh, met_thresh, vir_thresh)
    knnFunctions(arg_thresh, mge_thresh, met_thresh, vir_thresh)

def hcFunctions(arg_thresh, mge_thresh, met_thresh, vir_thresh): 
    """
    Runs hierarchical clustering for the main functional groups 
    Thresholded matrices used 
    Produces two plots: silhouette scores from HC and weighted Jaccard from KNN 
    """

    silhouette_arg = []
    silhouette_mge = []
    silhouette_met = []
    silhouette_vir = []

    silhouette_argT = []
    silhouette_mgeT = []
    silhouette_metT = []
    silhouette_virT = []


    ar_1 = []
    ar_2 = []
    ar_3 = []
    ar_4 = []

    maxL = 100

    runThrough = [arg_thresh, mge_thresh, met_thresh, vir_thresh]
    silThrough = [silhouette_argT, silhouette_mgeT, silhouette_metT, silhouette_virT]
    arThrough = [ar_1, ar_2, ar_3, ar_4]

    for i in range(2,maxL):
        
        for j in range(len(runThrough)):
            
            clust_temp, Z_temp = hc(runThrough[j], spec, "temp.nwk", i)
            c_temp = [k[0] for k in clust_temp]
        
            silThrough[j].append(silhouette_score(runThrough[j], c_temp))
            
            clust_ani, Z_ani = hc(a, spec, "ani_%i.nwk" %(i),i)
            c3 = [j[0] for j in clust_ani]
            
            arThrough[j].append(rand_score(c_temp, c3))

    hcFuncPlot(silhouette_argT,silhouette_mgeT,silhouette_metT,silhouette_virT)

def hcFuncPlot(silhouette_argT,silhouette_mgeT,silhouette_metT,silhouette_virT): 
    """
    Functions HC plot of silhouette scores
    """
    x = [i for i in range(2,maxL)]
    f,axes = plt.subplots(1,1, figsize=(8,6))
    plt.rcParams.update({'font.size': 30})

    sns.set_theme('paper')
    sns.set_style("white")

    sns.lineplot(x=x, y=silhouette_argT, label="ARG", linewidth=2, color="#2E86AB")
    axes.fill_between(x, silhouette_metT, silhouette_argT, alpha=.3, color="#2E86AB")

    sns.lineplot(x=x, y=silhouette_metT, label="MET", linewidth=2, color="#9D0D0D")
    axes.fill_between(x, silhouette_virT, silhouette_metT, alpha=.3, color="#9D0D0D")

    sns.lineplot(x=x, y=silhouette_virT, label="VIR", linewidth=2, color="#F15025")
    axes.fill_between(x, silhouette_mgeT, silhouette_virT, alpha=.3, color="#F15025")

    sns.lineplot(x=x, y=silhouette_mgeT, label="MGE", linewidth=2, color="#478E23")
    axes.fill_between(x, silhouette_syn, silhouette_mgeT, alpha=.3, color="#478E23")

    sns.lineplot(x=x, y=silhouette_syn, label="SYN", linewidth=2, color="#34003E")
    axes.fill_between(x, 0.25, silhouette_syn, alpha=.3, color="#34003E")
    sns.despine()

    plt.savefig("silhouetteScores_functional.png", dpi=700, bbox_inches="tight")
    print("Silhouette Scores: Functional")
    print("Mean, STD, Median")
    print("%0.3f, %0.3f, %0.3f" %(stats(silhouette_argT)))
    print("%0.3f, %0.3f, %0.3f" %(stats(silhouette_metT)))
    print("%0.3f, %0.3f, %0.3f" %(stats(silhouette_virT)))
    print("%0.3f, %0.3f, %0.3f" %(stats(silhouette_mgeT)))
    print()

def knnFunctions(arg_thresh, mge_thresh, met_thresh, vir_thresh):
    """
    Runs the KNN graphs for the functional thresholded groups 
    """
    runThrough = [arg_thresh, mge_thresh, met_thresh, vir_thresh]

    modularity_arg = []
    modularity_mge = []
    modularity_met = []
    modularity_vir = []

    maxL = 25

    wj_1 = []
    wj_2 = []
    wj_3 = []
    wj_4 = []

    dc_1 = []
    dc_2 = []
    dc_3 = []
    dc_4 = []

    runThrough = [arg_scaled, mge_scaled,
                met_scaled, vir_scaled]
    modThrough = [modularity_arg, modularity_mge,
                modularity_met, modularity_vir]
    wjThrough = [wj_1, wj_2, wj_3, wj_4]
    dcThrough = [dc_1, dc_2, dc_3, dc_4]

    for i in range(1,maxL):
        
        for j in range(len(runThrough)):
            
            A1 = kneighbors_graph(runThrough[j], i, mode='distance')
            G1 = nx.from_numpy_array(MinMaxScaler().fit_transform(A1.toarray()))
            c1 = nx.community.greedy_modularity_communities(G1)
            
            modThrough[j].append(nx.community.modularity(G1, c1))
            wjThrough[j].append(weightedJac(G1, G3))
            dcThrough[j].append(dice(G1, G3))
    
    print("Modularity: Functional")
    print("%0.3f, %0.3f, %0.3f" %(stats(modularity_arg)))
    print("%0.3f, %0.3f, %0.3f" %(stats(modularity_met)))
    print("%0.3f, %0.3f, %0.3f" %(stats(modularity_vir)))
    print("%0.3f, %0.3f, %0.3f" %(stats(modularity_mge)))
    print()

    print("Weighted Jaccard: Functional")
    print("Mean, STD, Median")
    print("%0.3f, %0.3f, %0.3f" %(stats(wj_1)))
    print("%0.3f, %0.3f, %0.3f" %(stats(wj_2)))
    print("%0.3f, %0.3f, %0.3f" %(stats(wj_3)))
    print("%0.3f, %0.3f, %0.3f" %(stats(wj_4)))
    print()

    wjFuncPlot(wj_1, wj_2, wj_3, wj_4)

def wjFuncPlot(wj_1, wj_2, wj_3, wj_4): 
    """
    Creates weighted jaccard plot for functions
    """
    x = [i for i in range(1,maxL)]

    f,axes = plt.subplots(1,1, figsize=(5,5))
    plt.rcParams.update({'font.size': 30})

    sns.set_theme('paper')
    sns.set_style("white")

    boxData = wj_1 + wj_2 + wj_3 + wj_4
    boxCols = ["ARG-ANI" for i in range(len(x))] + ["MGE-ANI" for i in range(len(x))] + ["MET-ANI" for i in range(len(x))] + ["VIR-ANI" for i in range(len(x))]
    data = pd.DataFrame()
    data['Weighted Jaccard'] = boxData
    data['Comparison'] = boxCols


    sns.stripplot(data=data, x="Comparison", y="Weighted Jaccard", hue=x+x+x+x, s=10)
    sns.despine()
    sns.move_legend(axes, "upper left", bbox_to_anchor=(1.5, 1))

    plt.savefig("weightedJaccard_functional.png", dpi=700, bbox_inches="tight")

def main(): 
    
    z,a,c = readInFiles() # READS IN FILES
    z_scaled, z_removed, z_thresh = scaling(z, a, c) # FORMS 3 SYNTENY METRICS

    hcOverall(z, a, c, spec, z_scaled, z_thresh, z_removed) # RUNS HC
    KNNGraphs(z, a, c, spec, z_scaled, z_thresh, z_removed) # RUNS KNN

    functions() # RUNS FUNCTIONS

    # TO CREATE FIGURE 3A/B OR FIGURE 4A/B
    # 3A/B: Produced by saving 


main()
