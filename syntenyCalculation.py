"""
cross counter
"""
import pandas as pd 
import numpy as np
import sys, os
from itertools import combinations
from sklearn.metrics.pairwise import cosine_similarity

def readGenomeData(species):

    locDF = pd.read_csv("/gpfs/data/isarkar/vramanan/cds_location.csv")
    locDF = locDF.drop_duplicates(subset='Species',keep='first')

    cols = ['Species','Locus','LocusTag','Start','Stop','Strand']
    geneDF = pd.read_csv("/gpfs/data/isarkar/vramanan/CDSDatabase_sampled.csv", usecols=cols)
    geneDF = geneDF[geneDF['Species'].isin(species)]
    geneDF = geneDF[geneDF['Locus'].isin(locDF['Locus'].values)]
    return geneDF

def skelGraph(currSpec, geneDF):
    df = geneDF[geneDF['Species']==currSpec].sort_values('Start')
    specNodes = []
    for i in range(len(df.index)):
        temp = df.iloc[i]
        specNodes.append(temp['LocusTag'])
    return specNodes

def getIndOrder(l, nodes1, nodes2):
    try: 
        loc1 = nodes1.index(l[0])
    except:
        loc1 = nodes1.index(l[1])
    try:
        loc2 = nodes2.index(l[1])
    except:
        loc2 = nodes2.index(l[0])
    return loc1, loc2

def crossCounter(edges, nodes1, nodes2, name):
    
    tot = list(set(list(combinations(edges, 2))))
    cross = 0
    totalEdges = 0 
    cosCross = 0
    cos = 0 
    for t in tot:
        tempCross = 0
        tempCos = 0 
        tempCosCross = 0 
        loc1a, loc1b = getIndOrder(t[0], nodes1, nodes2)
        loc2a, loc2b = getIndOrder(t[1], nodes1, nodes2)
        cosSim = cosine_similarity([[loc1a, loc1b], [loc2a, loc2b]])[0][1]
        cosDis = 1 - cosSim
        
        """
        if loc1a > loc2a: 
            if loc1b < loc2b:
                tempCosCross = cosSim
            else:
                tempCosCross = cosSim
        elif loc1a < loc2a:
            if loc1b > loc2b:
                tempCosCross = cosSim
            else:
                tempCosCross = cosSim
        """
        if loc1a > loc2a: 
            if loc1b < loc2b:
                cross += 1
        elif loc1a < loc2a:
            if loc1b > loc2b:
                cross += 1
        totalEdges += 1
        cosCross += cosSim

    if totalEdges < 10:
        print(name)
    return cosCross/len(tot), cross/len(tot)

def rearrangeCC(edges, nodes1, nodes2, curr):
    # add rearrangement here 
    crosses = []
    crossAvgs = []
    for e in edges:
        #print(e)
        t0 = e[0]
        t1 = e[1]
        if t0 in nodes1: 
            i0 = nodes1.index(t0)
            i1 = nodes2.index(t1)
            newOrd1 = nodes1[i0:len(nodes1)] + nodes1[0:i0]
            newOrd2 = nodes2[i1:len(nodes2)] + nodes2[0:i1]
        else:
            i0 = nodes2.index(t0)
            i1 = nodes1.index(t1)
            newOrd1 = nodes2[i0:len(nodes2)] + nodes2[0:i0]
            newOrd2 = nodes1[i1:len(nodes1)] + nodes1[0:i1]

        cosCross, crossAvg = crossCounter(edges, newOrd1, newOrd2, curr)
        crosses.append(cosCross)
        crossAvgs.append(crossAvg)
    print(curr)
    #print(crosses)
    minIndex = np.argmin(crosses)
    return crosses[minIndex], crossAvgs[minIndex]
    

def buildSkeleton(geneDF, edgeDict):
    skeletons = {}
    crossDict = {}
    counter = 0
    with open("crossCounter_v3.csv", 'w') as f:
        f.write("Pair,CosSim,CrossAvg\n")
        for curr in edgeDict:
            sp = curr.split(",")
            currSpec1 = sp[0]
            currSpec2 = sp[1]
            nodes1 = skelGraph(currSpec1, geneDF)
            nodes2 = skelGraph(currSpec2, geneDF)
            edges = edgeDict[curr]
            minCross, minNum = rearrangeCC(edges, nodes1, nodes2, curr)
            #break
            #cross, cos, cosCross = crossCounter(edges, nodes1, nodes2, curr)
            #print(curr, cross, cos, cosCross)
            f.write("%s,%s,%s\n" %(curr, minCross, minNum))
            f.flush()
            counter += 1
            #if counter == 100:
                #break
        
    return crossDict

def readEdgeData(df): 
    edgeDict = {}
    for i in range(len(df.index)):
        curr = df.iloc[i]
        pair = curr['SortedPair']
        loc1 = curr['qLocusTag']
        loc2 = curr['sLocusTag']
        if pair in edgeDict: 
            edgeDict[pair].append((loc1, loc2))
        else:
            edgeDict[pair] =[(loc1, loc2)]
    
    return edgeDict

def main():
    pairDF = pd.read_csv("CCdata/filteredEdges.csv")
    pairDF[['Spec1', 'Spec2']] = pairDF['SortedPair'].str.split(",", expand=True)[[0,1]]
    species = list(set(list(pairDF['Spec1'].unique()) + list(pairDF['Spec2'].unique())))
    geneDF = readGenomeData(species)
    edgeDict = readEdgeData(pairDF)
    print("building skeletons", flush=True)
    skeletons = buildSkeleton(geneDF, edgeDict)

    # get Edges from allDiffs 
    # subsample each set of edges to 10 maximum 
    # run cross counter

main()