"""
deltacon.py
--------------------------

Deltacon measure for graph distance, after:

Koutra, Danai, Joshua T. Vogelstein, and Christos Faloutsos. 2013. “Deltacon: A
Principled Massive-Graph Similarity Function.” In Proceedings of the 2013 SIAM
International Conference on Data Mining, 162–70. Society for Industrial and
Applied Mathematics. https://doi.org/10.1137/1.9781611972832.18.

author: Stefan McCabe
email: stefanmccabe at gmail dot com
Submitted as part of the 2019 NetSI Collabathon.

"""
import networkx as nx
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
from skbio.tree import TreeNode
from skbio import io
import matplotlib.pyplot as plt

def hc(X,labels,name):
    #X = MinMaxScaler().fit_transform(X)
    Z = linkage(X, 'complete')
    clusters = cut_tree(Z, n_clusters=5)

    tree = TreeNode.from_linkage_matrix(Z, labels)
    io.write(tree, format='newick', into="../Trees/"+name)
    
    fig = plt.figure(figsize=(15, 10))
    dn = dendrogram(Z, orientation='right', labels=labels)
    plt.rcParams.update({'font.size': 40})
    plt.show()
    
    return clusters

def dice(g1, g2):
    num = 2*len(list(set(g1.edges()).intersection(set(g2.edges()))))
    return num / (len(g1.edges()) + len(g2.edges()))

def weightedJac(G1, G):
    G1_dict = edgeIntoDict(G1.edges)
    G_dict = edgeIntoDict(G.edges)
    intersection = list(set(G1.edges()).intersection(set(G.edges())))
    num = 0.0
    denom = 0.0
    for l in intersection:
        num += min(G1_dict[l], G_dict[l])
        denom += max(G1_dict[l], G_dict[l])
    return num/denom


def dist(G1, G2, exact=True, g=None):
    """DeltaCon is based on the Matsusita between matrices created from fast
    belief propagation (FBP) on graphs G1 and G2.

    Because the FBP algorithm requires a costly matrix inversion, there
    is a faster, roughly linear, algorithm that gives approximate
    results.

    Parameters
    ----------

    G1, G2 (nx.Graph)
        two networkx graphs to be compared.

    exact (bool)
        if True, use the slower but exact algorithm (DeltaCon_0)

    g (int)
        the number of groups to use in the efficient algorithm. If
        exact is set to False but g is not set, the efficient algorithm
        will still behave like the exact algorithm, since each node is
        put in its own group.

    Returns
    -------

    dist (float)
        the distance between G1 and G2.

    References
    ----------

    .. [1] Koutra, Danai, Joshua T. Vogelstein, and Christos
           Faloutsos. 2013. "Deltacon: A Principled Massive-Graph
           Similarity Function." In Proceedings of the 2013 SIAM
           International Conference on Data Mining, 162–70. Society for
           Industrial and Applied
           Mathematics. https://doi.org/10.1137/1.9781611972832.18.

    """
    assert G1.number_of_nodes() == G2.number_of_nodes()
    N = G1.number_of_nodes()

    if not exact and g is None:
        g = N

    A1 = nx.to_numpy_array(G1)
    L1 = nx.laplacian_matrix(G1).toarray()
    D1 = L1 + A1

    A2 = nx.to_numpy_array(G2)
    L2 = nx.laplacian_matrix(G2).toarray()
    D2 = L2 + A2

    eps_1 = 1 / (1 + np.max(D1))
    eps_2 = 1 / (1 + np.max(D2))

    if exact:
        S1 = np.linalg.inv(np.eye(N) + (eps_1 ** 2) * D1 - eps_1 * A1)
        S2 = np.linalg.inv(np.eye(N) + (eps_2 ** 2) * D2 - eps_2 * A2)
    else:
        raise NotImplementedError(
            "The efficient algorithm is not "
            "implemented. Please use the exact "
            "algorithm."
        )

    def matusita_dist(X, Y):
        r"""Return the Matusita distance

        .. math::

            \sqrt{\sum_i \sum_j \left( \sqrt{X_{ij}} - \sqrt{Y_{ij}} \right)^{2}}


        between X and Y.
        """
        return np.sqrt(np.sum(np.square(np.sqrt(X) - np.sqrt(Y))))

    dist = matusita_dist(S1, S2)
    return 1/(1+dist)

def matusita_dist(X, Y):
    return np.sqrt(np.sum(np.square(np.sqrt(X) - np.sqrt(Y))))

def nodeattr(G1, G2):
    N = G1.number_of_nodes()

    A1 = nx.to_numpy_array(G1)
    L1 = nx.laplacian_matrix(G1).toarray()
    D1 = L1 + A1

    A2 = nx.to_numpy_array(G2)
    L2 = nx.laplacian_matrix(G2).toarray()
    D2 = L2 + A2

    eps_1 = 1 / (1 + np.max(D1))
    eps_2 = 1 / (1 + np.max(D2))

    S1 = np.linalg.inv(np.eye(N) + (eps_1 ** 2) * D1 - eps_1 * A1)
    S2 = np.linalg.inv(np.eye(N) + (eps_2 ** 2) * D2 - eps_2 * A2)

    W = []
    
    for i in range(len(A1)):
        diff = [a_i - b_i for a_i, b_i in zip(A1[i], A2[i])]
        summed = np.sum(np.abs(diff))
        if summed > 0:
            W.append(matusita_dist(S1[i], S2[i]))
            
    nodeAttrIndex = np.argsort(W)[::-1]
    
    E = []
    for v in range(len(nodeAttrIndex)):
        srcNode = nodeAttrIndex[v]
        r = [a_i - b_i for a_i, b_i in zip(A2[srcNode], A1[srcNode])]
        for k in range(len(r)):
            destNode = k
            edgeScore = W[srcNode] + W[destNode]
            if r[k] > 0:
                E.append([srcNode, destNode, '+', edgeScore])
            elif r[k] < 0:
                E.append([srcNode, destNode, '-', edgeScore])
    E.sort(key=lambda x: x[3], reverse=True)
                
    
    return nodeAttrIndex, W, E