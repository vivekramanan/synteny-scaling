"""
Post-Analysis
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
from sklearn.cluster import AffinityPropagation
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.spatial import procrustes
from skbio.stats.ordination import pcoa
import plotly.graph_objects as go
import plotly.express as px
from sklearn.metrics.cluster import contingency_matrix
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import silhouette_score
from sklearn.manifold import TSNE
import cmasher as cmr

sys.path.append("../synteny-scaling/")
import dendrogram_from_girvan_newman as dgn
import comparisons as comparisons
import networks as networks

def readingInData():
    sixteenSMashFile = ""
    wholeGenomeMashFile = ""
    coreFile = ""

    df = pd.read_csv(sixteenSMashFile,usecols=['Species1_Name','Species2_Name','MashDistance'])
    wgdf = pd.read_csv(wholeGenomeMashFile,usecols=['Genome1_Name','Genome2_Name','MashDistance'])
    wgdf.columns = ['Species1_Name','Species2_Name','MashDistance']

    coreDF = pd.read_csv(coreFile,names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])
    coreDF = coreDF[1:]

    spec = list(set(df['Species1_Name'].unique().tolist() + df['Species2_Name'].unique().tolist()))
    z = np.zeros((len(spec), len(spec)))

    for i in range(len(df.index)):
        temp = df.iloc[i]
        index1 = spec.index(temp['Species1_Name'])
        index2 = spec.index(temp['Species2_Name'])
        z[index1][index2] = temp['MashDistance']
        z[index2][index1] = temp['MashDistance']

    w = np.zeros((len(spec), len(spec)))

    for i in range(len(wgdf.index)):
        temp = wgdf.iloc[i]
        if temp['Species1_Name'] in spec and temp['Species2_Name'] in spec:
            index1 = spec.index(temp['Species1_Name'])
            index2 = spec.index(temp['Species2_Name'])
            w[index1][index2] = float(temp['MashDistance'])
            w[index2][index1] = float(temp['MashDistance'])   

    c = np.zeros((len(spec), len(spec)))
    for i in range(len(coreDF.index)):
        temp = coreDF.iloc[i]
        if temp['Species1_Name'] in spec and temp['Species2_Name'] in spec:
            index1 = spec.index(temp['Species1_Name'])
            index2 = spec.index(temp['Species2_Name'])
            if float(temp['CosSim']) >= 0.0: 
                c[index1][index2] = float(temp['CosSim'])
                c[index2][index1] = float(temp['CosSim'])
    c_rev = 1-c
    for i in range(len(c_rev)):
        c_rev[i][i] = 0.0 
    
    tax = pd.read_csv("../Data/species-class.csv")
    taxD = {}
    taxC = {}
    for i in range(len(tax.index)):
        if 'helicobacter' in tax.iloc[i]['species'] or 'campylobacter' in tax.iloc[i]['species']:
            taxD[tax.iloc[i]['species'].replace("'",'')] = 'Campylobacterota'
            taxC[tax.iloc[i]['species'].replace("'",'')] = 'Campylobacteria'
        else:
            taxD[tax.iloc[i]['species'].replace("'",'')] = tax.iloc[i]['phylum']
            taxC[tax.iloc[i]['species'].replace("'",'')] = tax.iloc[i]['class']
    
    return z,w,c_rev,taxD,taxC,spec

def scale(old, syn, length):
    scaled = np.zeros((length, length))
    for i in range(len(old)):
        for j in range(len(old)):
            if syn[i][j] != 0.0:
                scaled[i][j] = old[i][j] * syn[i][j]
    #scale = MinMaxScaler().fit_transform(scaled)
    return scaled

def hcAnalysis(clust_16, clust_16syn, clust_wg, z_scaled, wg_scaled, z, spec, taxD):
    c1 = [i[0] for i in clust_16]
    c2 = [i[0] for i in clust_16syn]
    c_wg = [i[0] for i in clust_wg]

    print("Sil 16S:", silhouette_score(z, c1))
    print("Sil 16SSyn:", silhouette_score(z_scaled, c2))
    print("Sil 16SWG:", silhouette_score(wg_scaled, c_wg))
    print()

    print("16s vs 16Syn:")
    print(adjusted_rand_score(c1, c2))
    print(contingency_matrix(c1, c2))
    print()

    print("16s vs 16WG:")
    print(adjusted_rand_score(c1, c_wg))
    print(contingency_matrix(c1, c_wg))
    print()

    # ADJUST FIT ACCORDINGLY BASED ON CONTINGENCY MATRIX 
    fit = {0:3,1:1,2:2,3:2,4:4}
    c3 = [fit[i] for i in c2]
    fit2 = {0:3,1:1,2:1,3:2,4:4}
    c4 = [fit2[i] for i in c_wg]

    phyla = {}
    for i in range(len(spec)):
        p = taxD[spec[i]]
        if p in phyla:
            phyla[p]+= 1
        else:
            phyla[p]=1
    
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
    
    PCAandTSNE(z, z_scaled, wg_scaled, c1, c3, c4, phylaColors, phylakey)

    return
    
def PCAandTSNE(z, z_syn, z_wg, c1, c3, c4, phylaColors, phylakey):
    # PCA and label clusters 
    f,axes = plt.subplots(1,2, figsize=(20,10))
    plt.rcParams.update({'font.size': 20})

    pca1 = PCA(n_components=2)
    pca1.fit(z)
    print(pca1.explained_variance_ratio_)

    pca2 = PCA(n_components=2)
    pca2.fit(z_syn)
    print(pca2.explained_variance_ratio_)

    for i in list(set(c1)):
        indices = [j for j, x in enumerate(c1) if x == i]
        x = pca1.components_[0][indices]
        y = pca1.components_[1][indices]
        axes[0].scatter(x, y, color=phylaColors[i],label=phylakey[i])
        
        indices = [j for j, x in enumerate(c3) if x == i]
        x = pca2.components_[0][indices]
        y = pca2.components_[1][indices]
        axes[1].scatter(x, y, color=phylaColors[i],label=phylakey[i])

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    #axes[0].set_title("16S Hierarchical Clusters")
    #axes[1].set_title("16s+Synteny Hierarchical Clusters")

    #axes[0].legend()
    #axes[1].legend()
    plt.savefig("../Figures/PCA.png", dpi=700, bbox_inches="tight")
    
    # TSNE
    f,axes = plt.subplots(1,3, figsize=(30,10))
    plt.rcParams.update({'font.size': 20})

    tsne1 = TSNE(n_components=2, random_state=0)
    proj1 = tsne1.fit_transform(z)

    tsne2 = TSNE(n_components=2, random_state=0)
    proj2 = tsne2.fit_transform(z_syn)

    tsne3 = TSNE(n_components=2, random_state=0)
    proj3 = tsne2.fit_transform(z_wg)

    for i in list(set(c1)):
        indices = [j for j, x in enumerate(c1) if x == i]
        x = proj1[:,0][indices]
        y = proj1[:,1][indices]
        axes[0].scatter(x, y, color=phylaColors[i],label=phylakey[i])
        
        indices = [j for j, x in enumerate(c3) if x == i]
        x = proj2[:,0][indices]
        y = proj2[:,1][indices]
        axes[1].scatter(x, y, color=phylaColors[i],label=phylakey[i])
        
        indices = [j for j, x in enumerate(c4) if x == i]
        x = proj3[:,0][indices]
        y = proj3[:,1][indices]
        axes[2].scatter(x, y, color=phylaColors[i],label=phylakey[i])

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        
    plt.savefig("../Figures/TSNE.png", dpi=700, bbox_inches="tight")
    plt.show()

def KNNgraphs(z, z_scaled, wg_scaled, taxC, taxD, spec):

    prism = px.colors.qualitative.Prism
    taxColor = {}
    s = list(set(taxD.values()))
    for i in range(len(s)):
        taxColor[s[i]] = prism[i]

    colors = cmr.take_cmap_colors('tab20', 20, return_fmt='hex')

    taxClassColor = {}
    c = list(set(taxC.values()))
    for i in range(len(c)):
        taxClassColor[c[i]] = colors[i]

    # KNN Graph: 16S Original (Phylum and Class)
    A1 = kneighbors_graph(z, 10, mode='distance')
    G1 = nx.from_numpy_array(MinMaxScaler().fit_transform(A1.toarray()))
    k1 = nx.kamada_kawai_layout(G1)
    networks.plotNetwork(G1, k1, taxColor, taxD, spec, "KNN Graph: 16S (Phylum)")
    networks.plotNetwork(G1, k1, taxClassColor, taxC, spec, "KNN Graph: 16S (Class)")

    # KNN Graph: 16S + Synteny (Phylum and Class)
    A = kneighbors_graph(z_scaled, 10, mode='distance')
    G = nx.from_numpy_array(MinMaxScaler().fit_transform(A.toarray()))
    k = nx.kamada_kawai_layout(G)
    networks.plotNetwork(G, k, taxColor, taxD, spec, "KNN Graph: 16S + Synteny (Phylum)")
    networks.plotNetwork(G, k, taxClassColor, taxC, spec, "KNN Graph: 16S + Synteny (Class)")

    # KNN Graph: 16S + Whole Genome (Phylum and Class)
    A2 = kneighbors_graph(wg_scaled, 10, mode='distance')
    G2 = nx.from_numpy_array(MinMaxScaler().fit_transform(A2.toarray()))
    k2 = nx.kamada_kawai_layout(G2)
    networks.plotNetwork(G2, k2, taxColor, taxD, spec, "KNN Graph: 16S+WG (Phylum)")
    networks.plotNetwork(G2, k2, taxClassColor, taxC, spec, "KNN Graph: 16S+WG (Phylum)")

    # Analysis: 16S vs 16S + Synteny
    G1_dict = edgeIntoDict(G1.edges)
    G_dict = edgeIntoDict(G.edges)
    intersect = len(set(G1.edges()).intersection(set(G.edges())))
    union = len(set(G1.edges()).union(set(G.edges())))
    print("Jaccard Similarity: ", intersect/union)
    print("Weighted Jaccard Similarity ", comparisons.weightedJac(G1, G))
    print("Dice Coefficient: ", comparisons.dice(G1, G))
    print("DeltaCon 16 vs Synt: ", 1-comparisons.dist(G, G1, exact=True, g=None))

    intersect = len(set(G1.edges()).intersection(set(G2.edges())))
    union = len(set(G1.edges()).union(set(G2.edges())))
    print("Jaccard Similarity 16 vs WG: ", intersect/union)
    print("Weighted Jaccard Similarity 16 vs WG: ", comparisons.weightedJac(G1, G2))
    print("Dice Coefficient 16 vs WG: ", comparisons.dice(G1, G2))
    print("DeltaCon 16 vs WG: ", 1-comparisons.dist(G2, G1, exact=True, g=None))

    # Girvan Newman Community Detection
    partitions = girvannewmanClust(G)
    bestGN = partitions[13] # CHOOSE PARTITIONS ACCORDING TO BEST ONE IN GIRVAN NEWMAN (can change)
    HCGNComparison(bestGN, z_scaled)

def HCGNComparison(bestGN, z_scaled):
    Z = linkage(z_scaled, 'complete')
    clusters = cut_tree(Z, n_clusters=15)

    hc_clust = [i[0] for i in clusters]
    gn_clust = []
    for i in range(len(hc_clust)):
        for j in range(len(bestGN)):
            if i in bestGN[j]:
                gn_clust.append(j)

    print(adjusted_rand_score(hc_clust, gn_clust))
    sns.heatmap(contingency_matrix(hc_clust, gn_clust))
    plt.savefig("../Figures/HCGNComparison.png", dpi=700, bbox_inches="tight")

def dcNodeEdgeAttr(G, G1, spec, taxC):
    # DeltaCon: Node and Edge Attr
    nodeAttrIndex, impact, sorted_E = comparisons.nodeattr(G, G1)
    for i in range(len(nodeAttrIndex)):
        print(spec[nodeAttrIndex[i]], impact[nodeAttrIndex[i]])
    for i in range(len(sorted_E)):
        print("%s --> %s: \t %s \t %s" %(spec[sorted_E[i][0]], spec[sorted_E[i][1]], sorted_E[i][2], sorted_E[i][3]))
    
    classNodeAttr = {}
    for i in range(len(nodeAttrIndex)):
        cl = taxC[spec[nodeAttrIndex[i]]]
        if cl in classNodeAttr:
            classNodeAttr[cl] += impact[nodeAttrIndex[i]]
        else:
            classNodeAttr[cl] = impact[nodeAttrIndex[i]]
    
    x = []
    y = []
    for w in sorted(classNodeAttr, key=classNodeAttr.get, reverse=True):
        x.append(w)
        y.append(classNodeAttr[w])

    f,axes = plt.subplots(1,1, figsize=(20,5))
    sns.barplot(x=x, y=y, color='tomato')
    plt.xticks(ticks=[],labels=[])
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['left'].set_visible(False)
    plt.savefig("../Figures/nodeAttr.png", dpi=700, bbox_inches="tight")
    plt.show()

    ea = np.zeros((len(x), len(x)))
    for i in range(len(sorted_E)):
        cl1 = taxC[spec[sorted_E[i][0]]]
        cl2 = taxC[spec[sorted_E[i][1]]]
        ea[x.index(cl1)][x.index(cl2)] += sorted_E[i][3]
        ea[x.index(cl2)][x.index(cl1)] += sorted_E[i][3]
    
    f,axes = plt.subplots(1,1, figsize=(20,10))
    sns.heatmap(MinMaxScaler().fit_transform(ea), cmap="rocket_r", xticklabels=[], yticklabels=[])
    plt.savefig("../Figures/edgeAttr.png", dpi=700, bbox_inches="tight")
    plt.show()

def edgeIntoDict(G_edges):
    d = {}
    edges = list(G_edges(data=True))
    for i in edges:
        d[(i[0], i[1])] = i[2]['weight']
    return d

def girvannewmanClust(G):
    # girvan newman dendrogram 
    nn = nx.number_of_nodes(G)
    ne = nx.number_of_edges(G)
    partitions = dgn.girvan_newman_partitions(G)

    # print the list of partitions
    for i, part in enumerate(partitions):
        print ('partition number %d:\t'%i, part)
        
    print ('\nN.B: number of partitions: %d\t is equal to\t (number of nodes - 1): %d - 1'
        %(len(partitions), nn))

    agglomerative_mat = dgn.agglomerative_matrix(G, partitions)
    print ('Agglomerative matrix:\n', agglomerative_mat)

    dendro = dendrogram(agglomerative_mat)

    bp_G, idx_bp_G = dgn.girvan_newman_best_partition(G, partitions)

    print ('The best community partition of the graph G is:\n\t', bp_G)
    print ('\nWhich is the one with index = %d among the partitions in the list "partitions".'
        %(idx_bp_G))

    # How many communities in the best partition?
    n_communities_bp = len(bp_G)
    print ('The best partition splits the graph G into %d communities.'%n_communities_bp)

    # Distance of the best partition from the ground level
    dist_bp = dgn.distance_of_partition(agglomerative_mat, n_communities_bp)
    print ('The best partition, which splits the graph G into %d communities, has a distance from the ground level equal to:\t %d'
        %(n_communities_bp, dist_bp))

    # Plot the dendrogram highlighting the communities which compose the best partition
    dendro_bp = dendrogram(agglomerative_mat, color_threshold=dist_bp)

    tree = TreeNode.from_linkage_matrix(agglomerative_mat, spec)
    io.write(tree, format='newick', into="../bacteria_genome/Trees/16SYN_girvannewman_dist.nwk")
    return partitions

def main():
    z,w,c_rev,taxD,taxC,spec = readingInData()

    ## HIERARCHICAL CLUSTERING

    # 16S original 
    clust_16 = comparisons.hc(z, spec, "16s_original.nwk")
    # 16s + Synteny 
    z_scaled = MinMaxScaler().fit_transform(np.dot(z, np.cov(c_rev)))
    clust_16syn = comparisons.hc(MinMaxScaler().fit_transform(np.dot(z, np.cov(c_rev))), spec, "DOT_16syn.nwk") 
    # WG + Synteny
    wg_scaled = wg_scaled = MinMaxScaler().fit_transform(np.dot(z, np.cov(w)))
    clust_wg = comparisons.hc(MinMaxScaler().fit_transform(np.dot(z, np.cov(w))), spec, "DOT_wgsyn.nwk")

    # HC Analysis
    hcAnalysis(clust_16, clust_16syn, clust_wg, z_scaled, wg_scaled, z, spec, taxD)

    # KNN Graphs
    KNNgraphs(z, z_scaled, wg_scaled, taxC, taxD, spec)
    
main()