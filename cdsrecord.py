"""
CDS Record Object Definition
Used in BLAST Formatting Script

Vivek Ramanan
"""

class CDSRecord():
    def __init__(self, start, stop, strand, sequence):
        self.start = start
        self.stop = stop
        self.strand = strand
        self.seq = sequence
        self.gene = ""
        self.description = ""
        self.locusTag = ""
        self.dbxref = ""
        self.proteinID = ""
        self.translation = ""
        self.locus = ""
        self.species = ""
