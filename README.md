# synteny-scaling
 
 Code and supplementary data for synteny scaling method presented in WIP manuscript.

 # 1 Ortholog Construction

 Constructing orthologs with command line BLAST, followed by filtering.
 db.fasta: database of sequences to be BLAST'ed against, which consisted of all genes chosen
 query.fasta: functional genes chosen from external database cohorts

 Example BLAST command:
     '''
     blastn -db db.fasta -query query.fasta -outfmt 10 -out queryout.csv
     '''
 Outformat 10 outputs a CSV file, which has been used for blastFormatting.py.

 Files:
 - blastFormatting.py: formats the BLAST CSV output files
 - cdsrecord.py: object of record

 Outputs of blastFormatting.py can be used however - I used them for step 2 after filtering for chosen species. 

 # 2 Synteny Measure Calculation

 Files:
 - syntenyCalculation.py

 Requires a set of edges and genome gene data (same one used previously in step 1, as a CSV however)
 Outputs a CSV with the synteny measure per species pair

 # 3 Post-Calculation Analysis

 Compounds post-calculation analysis of hierarchical clustering, KNN graphs, and other statistical techniques.

 Files:
 - postAnalysis_1.py: python script version synteny metric analysis 
        (includes creation of hierarchical clustering and KNN graphs figures)
 - postAnalysis_2.py: revised version of hierarchical clustering and KNN that produces plots that
        compare cluster cutoff and clustering quality, across other variables
 - taxonomy.py: calculates taxonomic data for pairs (requires taxDict.txt)
 - functions.py: runs the functions for postAnalysis.py
 - networks.py: network functions for postAnalysis.py
 - dendrograms_from_girvan_newman.py: written by Francesco Bonacina to find dendrogram after GN analysis
