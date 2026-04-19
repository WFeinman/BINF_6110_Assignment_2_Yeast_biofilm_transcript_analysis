# BINF_6110_Assignment_2_Yeast_biofilm_transcript_analysis

##Overal Goal:
Quantification and functional analysis of S. cerevisiae biofilm formation across development stages.

## Shell Script Execution:
Can be bypassed with the quant_output files contained in the repository, provided the .gtf reference file is downlowaded with:
-     wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz

Ensure the SRA toolkit and the Salmon package are properly installed and active. Last confirmed functional versions: salmon 1.10.2, SRA toolkit 3.4.1.

- Place the “Assignment_2_pseudoalignment_code.sh” script in desired file location. 
- Give the script execution privileges with:
-     chmod +x Assignment_2_pseudoalignment_code.sh

Execute the primary script with:
-     bash Assignment_2_pseudoalignment_code.sh

Once completed, a quant_output directory should be created in the current working directory. 


## R Script Execution:
Requires "s_cerevisae_biofilm_accession_lookup.csv" and the quant_output_directory in the current working directory

Open code in R-Studio and run all. Last confirmed functional with: R-4.5.1. BiocManager 1.30.27, ggplot2 4.0.2, tidyverse 2.0.0, AnnotationDBI 1.70.0, apeglm 1.3.0, clusterProfiler 4.16.0, DOSE 4.2.0, enrichplot 1.28.4, GenomicFeatures 1.60.0, pheatmap 1.0.13, txdbmaker 1.4.2, org.Sc.sgd.db 3.21.0. 
