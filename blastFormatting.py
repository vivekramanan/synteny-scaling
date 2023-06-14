"""
Reads in BLAST CSV, formats and filters based on filter thresholds
Also calculates intergenic distance between pairs of orthologs (not used in analysis)

Vivek Ramanan
"""
import pandas as pd
import numpy as np
import sys,os
from cdsrecord import CDSRecord
from scipy import stats
import taxonomy

def getDifference(blast, matchedData, filename):
    """
    Calculates the location difference between the first and second gene orthologs
    @param blast: dataframe 
    @param matchedData: 
    @param filename: 
    @return diff: 
    """
    diffs = []
    with open(filename, 'w') as f:
        f.write("qid,sid,difference\n")
        for i in range(len(blast.index)):
            temp = blast.iloc[i]
            qid = temp['qseqid'].split("|")[0]
            sid = temp['sseqid'].split("|")[0]
            if qid != sid: 
                diff = np.abs(matchedData.loc[qid]['Start'] - matchedData.loc[sid]['Start'])
                diffs.append(diff)
                f.write("%s,%s,%s\n" %(qid, sid, diff))
    print(stats.describe(diffs))
    return diff

def getMinTax(temp, taxonomyDict):
    """
    Finding the minimum shared taxonomic level between two species
    @param temp: dictionary that holds the first species (qSpecies) and second (sSpecies)
    @param taxonomyDict: dictionary with taxonomic information, organized with 
        keys of taxonomic level and values of the taxonomy name
    """
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
    """
    Finds the new start location if needed to be reversed
    @param start: starting location
    @param specLocus: particular locus of the species 
    @param locus: dictionary with all loci
    @return revStart: int of new start location 
    """
    end = locus[specLocus]['Stop']
    revStart = end - start 
    return revStart

def changeDirection(qid, sid, match, locus):
    """
    Checks if species genome data needs to be reverse (i.e. one + and one - strand or vice versa)
    @param qid: first species ID
    @param sid: second species ID
    @param match: dictionary with all the IDs and start and stop info
    @param locus: particular locus of the ID 
    @return new start of both qid and sid 
    """
    qStart = match[qid]['Start']
    sStart = match[sid]['Start']

    if match[qid]['Strand'] == "-" and match[sid]['Strand'] == "-":
        return qStart, sStart

    if match[qid]['Strand'] == "-":
        qStart = reverseStrand(qStart, match[qid]['Locus'], locus)
    
    if match[sid]['Strand'] == "-":
        sStart = reverseStrand(sStart, match[sid]['Locus'], locus)

    return qStart, sStart

def checkLocType(qid, sid, match, locus):
    """
    Checking whether the locus is genome or plasmid
    @param qid: first species ID
    @param sid: second species ID
    @param match: dictionary with all the IDs and start and stop info
    @param locus: particular locus of the ID 
    @return qGen and sGen: str with either "Genome" or "Plasmid" values
    """
    qGen = locus[match[qid]['Locus']]['GenType']
    sGen = locus[match[sid]['Locus']]['GenType']
    return qGen, sGen

def getDiff(blast, matchedData, f, taxDict, locus, sampled):
    """
    Writing the file with the blast information and intergenic distance (if needed)
    @param blast: dictionary with all the row information
    @param matchedData: data from the original entire gene file 
    @param f: file to write to 
    @param taxDict: taxonomy dictionary to calculate distance
    @param locus: dataframe of all genome loci (excluding plasmids)
    @param sampled: chosen GCFs to analyze
    """
    blast['qLocusTag'] = blast['qseqid'].str.split("|", expand=True)[0]
    blast['sLocusTag'] = blast['sseqid'].str.split("|", expand=True)[0]
    blast['qSpecies'] = blast['qseqid'].str.split("|", expand=True)[3]
    blast['sSpecies'] = blast['sseqid'].str.split("|", expand=True)[3]

    print(blast.shape)
    blast = blast[(blast['qSpecies'].isin(sampled['Species'])) & (blast['sSpecies'].isin(sampled['Species']))]
    print(blast.shape)

    tags = list(set(list(blast['qLocusTag'].unique()) + list(blast['sLocusTag'].unique())))
    matchTemp = matchedData[matchedData.index.isin(tags)]
    match = matchTemp.to_dict(orient='index')

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
                else:
                    diff = "NA"
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
    """
    Reading in a CSV which has information on which loci are genomes (first of the file)
    Var: location needs to be set accordingly
    """
    location = "cds_genloc.csv" # CHANGE ACCORDINGLY 
    locus = pd.read_csv(location, index_col="Locus").to_dict(orient="index")

    return locus

def underSampleBlast(f):
    """
    Removes any other GCFs that weren't chosen, and filters ONLY for 
        percent identical > 95  and length of alignment between 500 and 2500
    Writes to the same file as before with "sampled" suffix
    @param f: original file of written output
    """
    sampled = pd.read_csv("uniqueGCF_corr.csv")
    names = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart',
            'qend','sstart','send','evalue','bitscore']
    cols2 = ['qseqid','sseqid','length','pident']
    l = []
    for chunk in pd.read_csv(f, header=None, names=names, usecols=cols2, chunksize=10000):
        chunk['qSpecies'] = chunk['qseqid'].str.split("|", expand=True)[3]
        chunk['sSpecies'] = chunk['sseqid'].str.split("|", expand=True)[3]
        chunk = chunk[(chunk['qSpecies'].isin(sampled['Species'])) & (chunk['sSpecies'].isin(sampled['Species']))]
        chunk = chunk[(chunk['pident'] > 95.00) & (chunk['length'] <= 2500) & (chunk['length'] > 500)]
        l.append(chunk)
    df = pd.concat(l, ignore_index=True)

    df.to_csv(f[:-4]+"_sampled.csv", index=False)

    return 

def main():
    # use GCF file with all the organized GCF to species information and taxonomy info
    sampled = pd.read_csv("GCFFile.csv")
    taxDict = taxonomy.readTaxDict()

    # Reading in gene file with all the gene information 
    cols = ['LocusTag','Start','Stop','Locus','Strand']
    matchedData = pd.read_csv(sys.argv[2], usecols=cols, index_col="LocusTag")
    matchedData = matchedData.replace(r'\n',' ', regex=True) 
    matchedData = matchedData[~matchedData.index.duplicated(keep='first')]
    # print("read in CDS", flush=True)

    names = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart',
            'qend','sstart','send','evalue','bitscore']
    cols2 = ['qseqid','sseqid','length','pident']
    
    # reading in genomic loci information 
    locus = readInDB()

    # Opens output file for writing
    with open(sys.argv[3], 'w') as f:
        f.write("qLocusTag,sLocusTag,qSpecies,sSpecies,qStart,sStart,qGen,sGen,minTax,Distance\n")
        # Reads in BLAST file here 
        for chunk in pd.read_csv(sys.argv[1], chunksize=10000):
            getDiff(chunk, matchedData, f, taxDict, locus, sampled)
    
    # Samples chosen GCF, pident > 95, and alignment length filtering  
    underSampleBlast(sys.argv[3])

main()   