"""
Analysis of other functional cohorts: same style as postAnalysis.py just simplified

Files (all CSVs):
1. 16s Mash File
2. core gene cross counter
3. antibiotic resistance cross counter
4. mobile genetic elements cross counter
5. metabolism cross counter
6. virulence factor cross counter

"""

import pandas as pd
import networkx as nx
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics.cluster import contingency_matrix
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import silhouette_score
from sklearn.neighbors import kneighbors_graph
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.append("../synteny-scaling/")
import dendrogram_from_girvan_newman as dgn
import comparisons as comparisons
import networks as networks

def funcAdj(spec, df):
    """
    Filters the dataframe for the species names of choice and organizes
    into a 2D numpy array 
    @param spec: list of species names (with underscores)
    @param df: dataframe pandas
    @return the 2d numpy array 
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

def functionKNN(f, G1, taxColor, taxD, spec, title):
    """
    Runs all the KNN graphs and comparisons for the various functional cohorts
    Compares against the original 16S data 
    @param f: functional dataframe of choice
    @param G1: 16S graph 
    @param taxColor: coloring the taxa per phylum
    @param taxD: phylum taxonomy dictionary
    @param spec: list of species
    @param title: functional cohort name
    @return 16s+functional cohort synteny scaled data (2d matrix)
    """
    G2 = nx.from_numpy_array(f)
    G2.remove_nodes_from(list(nx.isolates(G2)))
    k2 = nx.spring_layout(G2)
    networks.plotIsoNetwork(G2, k2, taxColor, taxD, spec, "KNN Graph (Phylum): ", title)

    f_scaled = MinMaxScaler().fit_transform(np.dot(z, np.cov(f)))
    f_knn = kneighbors_graph(f_scaled, 10, mode='distance')
    fG = nx.from_numpy_array(f_knn)
    fG.remove_nodes_from(list(nx.isolates(fG)))
    f_layout = nx.kamada_kawai_layout(fG)
    networks.plotNetwork(fG, f_layout, taxColor, taxD, spec, "KNN Graph (Phylum): 16S + ", title)

    intersect = len(set(fG.edges()).intersection(set(G1.edges())))
    union = len(set(fG.edges()).union(set(G1.edges())))
    print("Jaccard Similarity: ", intersect/union)
    print("Weighted Jaccard Similarity ", comparisons.weightedJac(fG, G1))
    print("DeltaCon: ", 1-comparisons.dist(G1, fG, exact=True, g=None))
    print("Dice Coef: ", comparisons.dice(G1, fG))

    return f_scaled

def networks(spec, z, arg, mge, met, vir):
    """
    Runs all the KNN graph functions
    @param spec: list of species
    @param z: 16S 2d array 
    @param arg: antibiotic resistance df
    @param mge: mobile genetic elements df
    @param met: metabolism df
    @param vir: virulence df
    @return each of the synteny scaled 16S data per cohort
    """
    tax = pd.read_csv("../Data/species-class.csv")
    taxD = {}
    taxC = {}
    for i in range(len(tax.index)):
        taxD[tax.iloc[i]['species'].replace("'",'')] = tax.iloc[i]['phylum']
        taxC[tax.iloc[i]['species'].replace("'",'')] = tax.iloc[i]['class']
    prism = px.colors.qualitative.Prism
    taxColor = {}
    s = list(set(taxD.values()))
    for i in range(len(s)):
        taxColor[s[i]] = prism[i]
        
    alphabet = px.colors.qualitative.Alphabet
    taxClassColor = {}
    c = list(set(taxC.values()))
    for i in range(len(c)):
        taxClassColor[c[i]] = alphabet[i]

    A1 = kneighbors_graph(z, 10, mode='distance')
    G1 = nx.from_numpy_array(MinMaxScaler().fit_transform(A1.toarray()))

    # Each functional analysis
    arg_scaled = functionKNN(arg, G1, taxColor, taxD, spec, "Antibiotic Resistance")
    mge_scaled = functionKNN(mge, G1, taxColor, taxD, spec, "Mobile Genetic Elements")
    met_scaled = functionKNN(met, G1, taxColor, taxD, spec, "Metabolism Genes")
    vir_scaled = functionKNN(vir, G1, taxColor, taxD, spec, "Virulence Factors")

    return arg_scaled,mge_scaled,met_scaled,vir_scaled

def hc(z,arg_scaled,mge_scaled,met_scaled,vir_scaled):
    """
    Runs hierarchical clustering for each of the scaled synteny data
    @param z: 16S data original
    @param <function>_scaled: scaled functional 2d matrix
    @return returns c_<function>: clusters for each cohort per HC 
    """
    # 16S
    clust_16 = comparisons.hc(z, spec, "16s_original.nwk")
    # ARG
    clust_arg = comparisons.hc(arg_scaled, spec, "arg.nwk")
    # MGE
    clust_mge = comparisons.hc(mge_scaled, spec, "mge.nwk")
    # MET
    clust_met = comparisons.hc(met_scaled, spec, "met.nwk")
    # VIR
    clust_vir = comparisons.hc(vir_scaled, spec, "vir.nwk")

    c1 = [i[0] for i in clust_16]
    c_arg = [i[0] for i in clust_arg]

    print("ARG")
    print("Adj Rand Score: ", adjusted_rand_score(c1, c_arg))
    print(contingency_matrix(c1, c_arg))
    print("Silhouette Score: ", silhouette_score(arg_scaled, c_arg))
    print()

    c_met = [i[0] for i in clust_met]

    print("MET")
    print("Adj Rand Score: ", adjusted_rand_score(c1, c_met))
    print(contingency_matrix(c1, c_met))
    print("Silhouette Score: ", silhouette_score(met_scaled, c_met))
    print()

    c_mge = [i[0] for i in clust_mge]

    print("MGE")
    print("Adj Rand Score: ", adjusted_rand_score(c1, c_mge))
    print(contingency_matrix(c1, c_mge))
    print("Silhouette Score: ", silhouette_score(mge_scaled, c_mge))
    print()

    c_vir = [i[0] for i in clust_vir]

    print("VIR")
    print("Adj Rand Score: ", adjusted_rand_score(c1, c_vir))
    print(contingency_matrix(c1, c_vir))
    print("Silhouette Score: ", silhouette_score(vir_scaled, c_vir))

    fit_arg = {0:3,1:1,2:2,3:2,4:4}
    fit_met = {0:1,1:1,2:3,3:2,4:2}
    fit_mge = {0:1,1:1,2:3,3:2,4:4}
    fit_vir = {0:3,1:1,2:1,3:2,4:2}

    f_arg = [fit_arg[i] for i in c_arg]
    f_met = [fit_met[i] for i in c_met]
    f_mge = [fit_mge[i] for i in c_mge]
    f_vir = [fit_vir[i] for i in c_vir]

    m_arg = 0
    m_met = 0
    m_mge = 0
    m_vir = 0
    for i in range(len(c1)):
        if c1[i] == f_arg[i]:
            m_arg += 1
        if c1[i] == f_met[i]:
            m_met += 1
        if c1[i] == f_mge[i]:
            m_mge += 1
        if c1[i] == f_vir[i]:
            m_vir += 1
            
    print("ARG",m_arg/len(c1))
    print("VIR",m_vir/len(c1))
    print("MGE",m_mge/len(c1))
    print("MET",m_met/len(c1))
    return c_arg,c_mge,c_met,c_vir

def plots(arg_scaled,mge_scaled,met_scaled,vir_scaled,c_arg,c_mge,c_met,c_vir):
    """
    Plot creation for all the scaled data + clusters
        1. PCA for each 
        2. Comparisons for KNN 
        3. Comparisons for HC
    @param <function>_scaled: scaled synteny 2d matrix 
    @param c_<function>: clusters per function
    """
    phylakey = {0: "bacillota/actino", 
            1: "bacteroidota", 
            2: "pseudomonadota", 
            3: "spirochaetes", 
            4: "verrucomicrobia"}
    phylaColors = {0:'#440154',
                1:'#3b528b',
                2:'#21918c',
                3:'#5ec962',
                4:'#fde725'}

    # PCA and label clusters 
    f,axes = plt.subplots(2,2, figsize=(20,15))
    plt.rcParams.update({'font.size': 20})

    pca1 = PCA(n_components=2)
    pca1.fit(arg_scaled)
    print(pca1.explained_variance_ratio_)

    pca2 = PCA(n_components=2)
    pca2.fit(mge_scaled)
    print(pca2.explained_variance_ratio_)

    pca3 = PCA(n_components=2)
    pca3.fit(met_scaled)
    print(pca3.explained_variance_ratio_)

    pca4 = PCA(n_components=2)
    pca4.fit(vir_scaled)
    print(pca4.explained_variance_ratio_)

    for i in list(set(c1)):
        indices = [j for j, x in enumerate(c_arg) if x == i]
        x = pca1.components_[0][indices]
        y = pca1.components_[1][indices]
        axes[0][0].scatter(x, y, color=phylaColors[i],label=phylakey[i])
        
        indices = [j for j, x in enumerate(c_mge) if x == i]
        x = pca2.components_[0][indices]
        y = pca2.components_[1][indices]
        axes[0][1].scatter(x, y, color=phylaColors[i],label=phylakey[i])
        
        indices = [j for j, x in enumerate(c_met) if x == i]
        x = pca3.components_[0][indices]
        y = pca3.components_[1][indices]
        axes[1][0].scatter(x, y, color=phylaColors[i],label=phylakey[i])
        
        indices = [j for j, x in enumerate(c_vir) if x == i]
        x = pca4.components_[0][indices]
        y = pca4.components_[1][indices]
        axes[1][1].scatter(x, y, color=phylaColors[i],label=phylakey[i])

        axes[0][0].set_title("ARG")
        axes[0][1].set_title("MGE")
        axes[1][0].set_title("MET")
        axes[1][1].set_title("VIR")

        axes[0][1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        #plt.savefig("../bacteria_genome/Figures/functionPCA.png", dpi=700, bbox_inches="tight")
        plt.show()
    
    # BAR GRAPH FOR KNN - hard coded, change accordingly
    x = ['Jaccard','Weighted\nJaccard','DICE','DeltaCon']
    y_cor = [0.665, 0.798, 0.799, 0.900]
    y_arg = [0.590, 0.499, 0.742, 0.926]
    y_met = [0.615, 0.569, 0.762, 0.922]
    y_vir = [0.629, 0.586, 0.772, 0.920]
    y_mge = [0.668, 0.596, 0.801, 0.921]

    y = y_arg + y_vir + y_mge + y_met + y_cor
    x = x*5
    z = []
    vals = ['ARG','VIR','MGE','MET','COR']
    for v in vals:
        z += [v for i in range(4)]
    
    f,axes = plt.subplots(1,1, figsize=(20,5))
    plt.rcParams.update({'font.size': 20})
    df = pd.DataFrame()
    df['Measure'] = x
    df['Value'] = y
    df['Type'] = z
    w = 0.7
    sns.barplot(data=df, x="Type", y="Value", hue="Measure", palette=["#80DED9","#307351","#7F7CAF","#310A31"])
    sns.set(font_scale=1.5)
    sns.set_style("white")
    sns.set_context("paper")

    axes.legend(bbox_to_anchor=(1, 1.05), prop={'size': 30})
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    plt.savefig("../Figures/networkComparisons.png", dpi=700, bbox_inches="tight")
    plt.show()

    # BAR GRAPH FOR HC - hard coded, change accordingly
    x = ['Adjusted\nRand','Silhouette','Manual\nMatch']
    y_cor = [0.42, 0.40, 0.78]
    y_arg = [0.739, 0.461, 0.907]
    y_met = [0.3382, 0.428, 0.834]
    y_mge = [0.466, 0.453, 0.836]
    y_vir =  [0.409, 0.418, 0.825]

    y = y_arg + y_vir + y_mge + y_met + y_cor
    x = x*5
    z = []
    vals = ['ARG','VIR','MGE','MET','COR']
    for v in vals:
        z += [v for i in range(3)]

    f,axes = plt.subplots(1,1, figsize=(20,5))
    df = pd.DataFrame()
    df['Measure'] = x
    df['Value'] = y
    df['Type'] = z
    w = 0.7
    plt.rcParams.update({'font.size': 20})
    sns.barplot(data=df, x="Type", y="Value", hue="Measure", palette=["#912F40","#202C39","#7D938A"])
    sns.set(font_scale=1.5)
    sns.set_style("white")
    sns.set_context("paper")

    axes.legend(bbox_to_anchor=(1, 1.05), prop={'size': 30})
    #axes.legend(bbox_to_anchor=(1.2, 1))
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    plt.savefig("../Figures/hcComparisons.png", dpi=700, bbox_inches="tight")
    plt.show()

def boxplot(argdf, virdf, mgedf, metdf, coreDF):
    """
    Boxplot for the synteny data per cohort distribution
    @param <function>df: dataframe for each 
    """
    # BOX PLOT
    order = {'ARG':argdf, 'VIR': virdf, 'MGE': mgedf, 'MET': metdf, 'COR':coreDF}
    vals = []
    types = []
    for o in order: 
        tempDF = order[o]
        types += [o for i in range(len(tempDF))]
        vals += tempDF['CosSim'].values.tolist()
    boxDF = pd.DataFrame()
    boxDF['Type'] = types
    boxDF['Values'] = vals
    boxDF['Values'] = pd.to_numeric(boxDF['Values'])
    sns.set_context('paper')
    f,axes = plt.subplots(1,1, figsize=(20,5))
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    plt.rcParams.update({'font.size': 20})
    sns.boxplot(data=boxDF, x='Type', y='Values')
    plt.savefig("../Figures/funcBox.png", dpi=700, bbox_inches="tight")

def main():
    df = pd.read_csv("16S_MASH_FILE",usecols=['Species1_Name','Species2_Name','MashDistance'])
    spec = list(set(df['Species1_Name'].unique().tolist() + df['Species2_Name'].unique().tolist()))
    z = np.zeros((len(spec), len(spec)))
    for i in range(len(df.index)):
        temp = df.iloc[i]
        index1 = spec.index(temp['Species1_Name'])
        index2 = spec.index(temp['Species2_Name'])
        z[index1][index2] = temp['MashDistance']
        z[index2][index1] = temp['MashDistance']

    mgedf = pd.read_csv("<MGE_File>", names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])
    mgedf = mgedf[1:]
    mge = funcAdj(spec, mgedf)

    metdf = pd.read_csv("<MET_FILE>", names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])
    metdf = metdf[1:]
    met = funcAdj(spec, metdf)

    virdf = pd.read_csv("<VIR_FILE>", names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])
    virdf = virdf[1:]
    vir = funcAdj(spec, virdf)

    argdf = pd.read_csv("<ARG_FILE>", names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])
    argdf = argdf[1:]
    arg = funcAdj(spec, argdf)

    coredf = pd.read_csv("CORE_FILE", names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])

    # KNN 
    arg_scaled,mge_scaled,met_scaled,vir_scaled = networks(spec, z, arg, mge, met, vir)

    # HC
    c_arg,c_mge,c_met,c_vir = hc(z,arg_scaled,mge_scaled,met_scaled,vir_scaled)

    # PLOTS
    plots(arg_scaled,mge_scaled,met_scaled,vir_scaled,c_arg,c_mge,c_met,c_vir)
    boxplot(argdf, virdf, mgedf, metdf, coredf)

main()