"""
ccross counter of ARGS specifically
"""
import pandas as pd 
import sys, os
from itertools import combinations
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np

def readGenomeData(species):

    locDF = pd.read_csv("/gpfs/data/isarkar/vramanan/cds_location.csv")
    locDF = locDF.drop_duplicates(subset='Species',keep='first')

    cols = ['Species','Locus','LocusTag','Start','Stop','Strand']
    geneDF = pd.read_csv("/gpfs/data/isarkar/vramanan/CDSDatabase_unique.csv", usecols=cols)
    geneDF = geneDF[geneDF['Species'].isin(species)]
    #geneDF = geneDF[geneDF['Locus'].isin(locDF['Locus'].values)]
    geneDF['splitLocus'] = geneDF['LocusTag'].str.split("_", expand=True)[1]
    return geneDF

def skelGraph(currSpec, geneDF):
    df = geneDF[geneDF['Species']==currSpec].sort_values('Start')
    specNodes = []
    for i in range(len(df.index)):
        temp = df.iloc[i]
        specNodes.append(temp['splitLocus'])
    return specNodes

def getIndOrder(l, nodes1, nodes2):
    if l[0] in nodes1: 
        loc1 = nodes1.index(l[0])
    elif l[0] in nodes2:
        loc1 = nodes2.index(l[0])
    else: 
        loc1 = None

    if l[1] in nodes1: 
        loc2 = nodes1.index(l[1])
    elif l[1] in nodes2:
        loc2 = nodes2.index(l[1])
    else: 
        print(l[1])
        loc2 = None

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
        try:
            cosSim = cosine_similarity([[loc1a, loc1b], [loc2a, loc2b]])[0][1]
        except Exception as e:
            #print(e)
            #print([[loc1a, loc1b], [loc2a, loc2b]])
            continue
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

    #if totalEdges < 10:
        #print(name)
    return cosCross/len(tot), cross/len(tot)

def rearrangeCC(edges, nodes1, nodes2, curr):
    # add rearrangement here 
    crosses = []
    crossAvgs = [0]
    #print(nodes1)
    #print(nodes2)
    print(edges)
    for edge in edges:
        #print(e)
        t0 = edge[0]
        t1 = edge[1]
        if t0 in nodes1: 
            try:
                i0 = nodes1.index(t0)
                i1 = nodes2.index(t1)
                newOrd1 = nodes1[i0:len(nodes1)] + nodes1[0:i0]
                newOrd2 = nodes2[i1:len(nodes2)] + nodes2[0:i1]
            except Exception as e:
                print(1, e)
                continue
        else:
            try:
                i0 = nodes2.index(t0)
                i1 = nodes1.index(t1)
                newOrd1 = nodes2[i0:len(nodes2)] + nodes2[0:i0]
                newOrd2 = nodes1[i1:len(nodes1)] + nodes1[0:i1]
            except Exception as e:
                print(2, e)
                continue

        cosCross, crossAvg = crossCounter(edges, newOrd1, newOrd2, curr)
        if cosCross:
            crosses.append(cosCross)
        if crossAvg:
            crossAvgs.append(crossAvg)
    #print(curr)
    #print(crosses)
    #print(crossAvgs)
    #minIndex = np.argmin(crosses)
    if len(crosses)>0 and len(crossAvgs)>0:
        return np.median(crosses), np.max(crossAvgs)
    else:
        print(curr, len(crosses), len(crossAvgs))
        return 0,0

def buildSkeleton(geneDF, edgeDict, name):
    skeletons = {}
    crossDict = {}
    counter = 0
    with open("crossCounter_"+name+".csv", 'w') as f:
        f.write("Pair,CosSim,CrossAvg\n")
        for curr in edgeDict:
            if len(edgeDict[curr]) > 1:
                sp = curr.split(",")
                currSpec1 = sp[0]
                currSpec2 = sp[1]
                nodes1 = skelGraph(currSpec1, geneDF)
                nodes2 = skelGraph(currSpec2, geneDF)
                edges = edgeDict[curr]
                minCross, minNum = rearrangeCC(edges, nodes1, nodes2, curr)
                f.write("%s,%s,%s\n" %(curr, minCross, minNum))
                f.flush()
                counter += 1
                #break
                #if counter == 5:
                    #break
        
    return crossDict

def readEdgeData(df): 
    edgeDict = {}
    for i in range(len(df.index)):
        curr = df.iloc[i]
        pair = curr['Pair']
        #loc1 = curr['qLocusTag']
        #loc2 = curr['sLocusTag']
        loc1 = curr['qLocusTag'].split("_")[1]
        loc2 = curr['sLocusTag'].split("_")[1]
        if pair in edgeDict: 
            edgeDict[pair].append((loc1, loc2))
        else:
            edgeDict[pair] =[(loc1, loc2)]
    
    return edgeDict

def main():
    
    cols = ['qLocusTag','sLocusTag','qSpecies','sSpecies','qGen','sGen']

    name = sys.argv[2]
    filtDF = pd.read_csv(sys.argv[1])
    filtDF = filtDF[(filtDF['qGen']=="Genome") & (filtDF['sGen']=="Genome")]
    species = list(set(list(filtDF['qSpecies'].unique()) + list(filtDF['sSpecies'].unique())))
    geneDF = readGenomeData(species)
    edgeDict = readEdgeData(filtDF)
    print("building skeletons", flush=True)
    skeletons = buildSkeleton(geneDF, edgeDict, name)

    # get Edges from allDiffs 
    # subsample each set of edges to 10 maximum 
    # run cross counter

main()