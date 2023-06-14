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

def readingInData():
    df = pd.read_csv("../bacteria_genome/Data/16s_mash2.csv",usecols=['Species1_Name','Species2_Name','MashDistance'])
    wgdf = pd.read_csv("../bacteria_genome/Data/mashUnique.csv",usecols=['Genome1_Name','Genome2_Name','MashDistance'])
    wgdf.columns = ['Species1_Name','Species2_Name','MashDistance']

    coreDF = pd.read_csv("../bacteria_genome/Data/crossCounter_CORE4.csv",names=['Species1_Name','Species2_Name','CosSim','CrossAvg'])
    coreDF = coreDF[1:]

    spec = list(set(df['Species1_Name'].unique().tolist() + df['Species2_Name'].unique().tolist()))
    z = np.zeros((len(spec), len(spec)))

    for i in range(len(df.index)):
        temp = df.iloc[i]
        index1 = spec.index(temp['Species1_Name'])
        index2 = spec.index(temp['Species2_Name'])
        z[index1][index2] = temp['MashDistance']
        z[index2][index1] = temp['MashDistance']

    

def main():
    

main()