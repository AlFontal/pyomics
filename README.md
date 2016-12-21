Repository where the data and scripts for the Pyomics group in Advanced Bioinformatics will be stored. 

The python script for the automated pipeline of the tuxedo protocol tools can be found in:

python_pipeline_tools.py

Before the clustering analysis script could be used we needed to convert the data that was obtained from stringtie to make it suitable to analyze in R:

StringtieToCluster.R


The script for the clustering analysis can be found in:

Clustering_Rscript_BIS1.R

The script for the Differential expression analysis can be found in:

ballgown.R


There are 2 files with genes that were obtained from both analysis:

from the clustering analysis: Cluster_Interesting_genes.csv

from the DE analysis: ballgown_DE_genes.csv

From these files we only extracted the CRO_names, and were put in 2 new files:

DE_genes_ballgown.txt & DE_genes_clustering.txt

We then searched for the overlapping genes between the 2 and this is outputtet in:

overlapping_genes.csv

The functional annotation was linked to these genes and is found in :

overlapping_genes_with_annotation.txt


After these overlapping genes were obtained we build or own BLAST database by using the protein sequences in:

protseq.fsa

