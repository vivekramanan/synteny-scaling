"""
Read in BLAST CSV
"""
import pandas as pd
import numpy as np
import sys,os
from cdsrecord import CDSRecord
from scipy import stats
import taxonomy

def getDifference(blast, matchedData, filename):
    diffs = []
    with open(filename, 'w') as f:
        f.write("qid,sid,difference\n")
        for i in range(len(blast.index)):
            temp = blast.iloc[i]
            qid = temp['qseqid'].split("|")[0]
            sid = temp['sseqid'].split("|")[0]
            if qid != sid: 
                #diff = np.abs(matchedData[qid].start - matchedData[sid].start)
                diff = np.abs(matchedData.loc[qid]['Start'] - matchedData.loc[sid]['Start'])
                diffs.append(diff)
                f.write("%s,%s,%s\n" %(qid, sid, diff))
    print(stats.describe(diffs))
    return diff

def getMinTax(temp, taxonomyDict):
    qSpec = temp['qSpecies'].split("_")[0:2]
    sSpec = temp['sSpecies'].split("_")[0:2]
    if qSpec[0] == sSpec[0]:
        if qSpec[1] == sSpec[1]:
            return "Species"
        else:
            return "Genus"
    else:
        qData = taxonomyDict[qSpec[0].capitalize()]
        sData = taxonomyDict[sSpec[0].capitalize()]
        best = "Kingdom"
        if "phylum" in qData and "phylum" in sData:
            if qData['phylum'] == sData['phylum']:
                best = "Phylum"
        if "class" in qData and "class" in sData:
            if qData['class'] == sData['class']:
                best = "Class"
        if "order" in qData and "order" in sData:
            if qData['order'] == sData['order']:
                best = "Order"
        if "family" in qData and "family" in sData:
            if qData['family'] == sData['family']:
                best = "Family"
        return best
    return "None"

def reverseStrand(start, specLocus, locus):
    end = locus[specLocus]['Stop']
    #print(start, specLocus, end, flush=True)
    revStart = end - start 
    return revStart

def changeDirection(qid, sid, match, locus):
    qStart = match[qid]['Start']
    sStart = match[sid]['Start']
    #print(qStart, sStart, flush=True)
    #print(match[qid]['Strand'], match[sid]['Strand'], flush=True)

    if match[qid]['Strand'] == "-" and match[sid]['Strand'] == "-":
        return qStart, sStart

    if match[qid]['Strand'] == "-":
        qStart = reverseStrand(qStart, match[qid]['Locus'], locus)
    
    if match[sid]['Strand'] == "-":
        sStart = reverseStrand(sStart, match[sid]['Locus'], locus)
    
    #print(qStart, sStart, flush=True)
    #print()

    return qStart, sStart

def checkLocType(qid, sid, match, locus):
    qGen = locus[match[qid]['Locus']]['GenType']
    sGen = locus[match[sid]['Locus']]['GenType']
    return qGen, sGen

def relativePosition(qid, sid, match, locus):
    qRel = float(match[qid]['Start']) / float(locus[match[qid]['Locus']]['Stop'])
    sRel = float(match[sid]['Start']) / float(locus[match[sid]['Locus']]['Stop'])
    return qRel, sRel

def getDiff(blast, matchedData, f, taxDict, locus, sampled):
    blast['qLocusTag'] = blast['qseqid'].str.split("|", expand=True)[0]
    blast['sLocusTag'] = blast['sseqid'].str.split("|", expand=True)[0]
    blast['qSpecies'] = blast['qseqid'].str.split("|", expand=True)[3]
    blast['sSpecies'] = blast['sseqid'].str.split("|", expand=True)[3]

    print(blast.shape)
    blast = blast[(blast['qSpecies'].isin(sampled['Species'])) & (blast['sSpecies'].isin(sampled['Species']))]
    print(blast.shape)

    tags = list(set(list(blast['qLocusTag'].unique()) + list(blast['sLocusTag'].unique())))
    matchTemp = matchedData[matchedData.index.isin(tags)]
    #print(matchTemp.index)
    #matchTemp.index.drop_duplicates()
    #print(matchTemp.index)
    match = matchTemp.to_dict(orient='index')
    #print("ready to go", flush=True)

    for i in range(len(blast.index)):
        temp = blast.iloc[i]
        qid = temp['qLocusTag']
        sid = temp['sLocusTag']
        qgenus = temp['qSpecies'].split("_")[0]
        pident = float(temp['pident'])
        length = int(temp['length'])
        try: 
            sgenus = temp['sSpecies'].split("_")[0]
        except:
            print("Error:",temp)
            continue
        qGen, sGen = checkLocType(qid, sid, match, locus)
        genType = qGen + "-" + sGen
        if qid != sid: # given pident > 95 and length < 2500
            try: 
                qStart,sStart = changeDirection(qid, sid, match, locus)
                if qGen == sGen:
                    diff = np.abs(qStart - sStart)
                    #qRel, sRel = relativePosition(qid, sid, match, locus)
                    #relDiff = np.abs(qRel - sRel)
                else:
                    diff = "NA"
                    #relDiff = "NA"
            except:
                print(qid, sid, "fail")
                continue
            withinGenus = 0
            if qgenus == sgenus:
                withinGenus = 1
            minTax = getMinTax(temp, taxDict)
            f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(qid, sid, temp['qSpecies'], temp['sSpecies'], qStart, sStart, qGen, sGen, minTax, diff))
    return

def readInDB():
    location = "/gpfs/data/isarkar/vramanan/cds_genloc.csv"
    #strand = "/gpfs/data/isarkar/vramanan/cds_locus.csv"

    locus = pd.read_csv(location, index_col="Locus").to_dict(orient="index")
    #strand = pd.read_csv(strand, index_col="LocusTag").to_dict(orient="index")

    return locus

def main():
    # use sampled
    sampled = pd.read_csv("/users/vramanan/bacteria_genome_pipeline/sampledGCF_human.csv")

    taxDict = taxonomy.readTaxDict()
    cols = ['LocusTag','Start','Stop','Locus','Strand']
    matchedData = pd.read_csv(sys.argv[2], usecols=cols, index_col="LocusTag")
    matchedData = matchedData.replace(r'\n',' ', regex=True) 
    matchedData = matchedData[~matchedData.index.duplicated(keep='first')]
    #matchedData = matchedData.loc[matchedData.index.drop_duplicates()]

    names = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart',
            'qend','sstart','send','evalue','bitscore']
    cols2 = ['qseqid','sseqid','length','pident']
    print("read in CDS", flush=True)
    locus = readInDB()

    with open(sys.argv[3], 'w') as f:
        f.write("qLocusTag,sLocusTag,qSpecies,sSpecies,qStart,sStart,qGen,sGen,minTax,Distance\n")
        # header=None, names=names, usecols=cols2,
        for chunk in pd.read_csv(sys.argv[1], chunksize=10000):
            getDiff(chunk, matchedData, f, taxDict, locus, sampled)

main()   