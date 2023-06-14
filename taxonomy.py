"""
Taxonomy Reference
"""
import pandas as pd 
import sys

def readTaxDict():
    """
    Reading in taxonomy dictionary information 
        Some manual changes added based on NCBI taxonomy update (helicobacter, campylobacter)
        Current genera at the time of code still refers to old bacterial phyla 
    CHANGE input pathway accordingly 
    """
    taxonomyDict = {}
    with open("taxonomyDict.txt", 'r') as f:
        for line in f:
            split = line.rstrip().split("\t")
            key = split[0].split(" ")[0]
            vals = split[1:]
            keyDict = {}
            for value in vals:
                splitVal = value.split(":")
                keyDict[splitVal[0]] = splitVal[1]
            taxonomyDict[key] = keyDict
    taxonomyDict["Sarcina"] = {"phylum":"Firmicutes","class":"Clostridia","order":"Clostridiales","genus":"Sarcina","resolution":"genus"}
    taxonomyDict["Escherichia"] = {"phylum":"Proteobacteria","class":"Gammaproteobacteria","order":"Enterobacterales","family":"Enterobacteriaceae","genus":"Escherichia","resolution":"genus"}
    taxonomyDict['Campylobacter'] = {"phylum":"Campylobacterota","class":"Campylobacteria","order":"Campylobacterales","family":"Campylobacteraceae","genus":"Campylobacter","resolution":"genus"}
    taxonomyDict['Helicobacter'] = {"phylum":"Campylobacterota","class":"Campylobacteria","order":"Campylobacterales","family":"Helicobacteraceae","genus":"Helicobacter","resolution":"genus"}
    print(len(taxonomyDict))
    return taxonomyDict

def read16s():
    """
    Reading in 16S MASH distance output CSV file 
    Change pathway accordingly 
    """
    path = "pair16s.csv"
    df = pd.read_csv(path)
    df['Pair'] = df['spec1'] + "-" + df['spec2']
    df = df.set_index('Pair')
    return df

def get16sDist(df, t1, t2):
    """
    Checks if there is a 16S distance information
        Based on MASH algorithm, distances above a certain threshold automatically go to 1 
    @param df: 16s data
    @param t1: species 1
    @param t2: species 2
    @return 16S distance if found, else None
    """
    lookfor = t1 + "-" + t2
    try:
        found = df.loc[lookfor]['dist']
        return found
    except:
        print(lookfor, "no 16s info")
        return None

def calculateTax(t1, t2):
    """
    Calculates the minimum shared taxonomic level between two species
    @param t1: a dictionary of all the taxonomic levels for the first species
    @param t2: dictionary of taxonomic levels for second species
    @return best: the minimum shared taxonomic level 
    """
    best = None
    try:
        if t1['family'] != t2['family']:
            if t1['order'] != t2['order']:
                if t1['class'] != t2['class']:
                    if t1['phylum'] != t2['phylum']:
                        best="kingdom"
                    else:
                        best="phylum"
                else:
                    best = "class"
            else:
                best = "order"
        else:
            best = "family"
        return best
    except:
        try:
            if t1['order'] != t2['order']:
                if t1['class'] != t2['class']:
                    if t1['phylum'] != t2['phylum']:
                        best="kingdom"
                    else:
                        best="phylum"
                else:
                    best = "class"
            else:
                best = "order"
            return best
        except:
            print(t1, t2)
    return best