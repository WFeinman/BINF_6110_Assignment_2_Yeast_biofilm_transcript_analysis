---
noteID: fc0b3f74-f7e8-401c-9afc-3eb634b8dc72
---
# Transcription Analysis of *S.cerevisiae* Biofilm Development

## Introduction

*Saccharomyces cerevisiae* is an ecologically important species of yeast for a wide variety of industries, alcohol production in particular. While fermentation is a fundamental part of this process, many of the nuances of wine are dependent on the metabolic transition of flor yeast over the course of biofilm formation. To briefly summarize this process, yeast in oxygen-deprived, sugar rich conditions require an oxygen-independent metabolism to survive, which produces alcohol as a byproduct. These yeast populations are a distributed single-cell structure, maximizing surface area for nutrient absorption. When fermentation transitions to wine aging, sugar sources have been exhausted and a controlled amount of air is deliberately introduced. Yeast turns to oxygen-dependent processes, managing nutrient scarcity and alcohol metabolism as the remaining food source. Yeast population structure shifts to form biofilms called velum; forming mats on the surface of the wine to ensure air access for the colony, exchange limited nutrients, and provide protection from environmental stressors. The ongoing metabolic byproducts from this process are important determinants of the final taste and quality of a given bottle of aged wine. (Fanning & Mitchell, 2012; Mardanov et al., 2020)

The particulars of that process, however, vary wildly between yeast strains, and even between wineries. This process is highly complex, and subject to many shifts in input from strain, time, population, and environmental factors. Gaining a greater biological understanding of why and how these processes occur, in addition to being an interesting biological-and-metabolic pathway question, would be greatly beneficial to ensuring more efficient wine production and good drinks. (Esteve-Zarzoso et al., 2001)

Given that these shifts are not solely dependent on species or yeast strain, a transcriptome approach will be more informative as a functional explanation than comparing genetics. As the metabolism changes over time, doing so requires finding how the transcriptional drivers of that process likewise change over time, especially during the transition period of biofilm formation. 

In order to do so, this paper will be studying transcriptome data from a prior study on yeast biofilm formation, attempting to find transcripts which show significant activity change between growth states, and identify the biological processes associated with them.


## Methods:

Transcriptome data is downloaded from the Mardanov et al 2020 study on *S. cerevisiae* strain I-329, an industrial flor yeast. The data consists of 9 total sets of transcriptome sequences, with 3 replicates for each landmark stage of biofilm development (early biofilm, thin biofilm, mature biofilm). The data is downloaded as fastq files using the SRA toolkit from NCBI. (Mardanov et al., 2020; NCBI, 2026)

FastQC inspection was performed on each fastq file to ensure read quality.

Once downloaded, the "Salmon" pseudoaligner shell software is used to generate transcript counts from each transcriptome, using the reference *S. cerevisiae* sequence from the National Institute of Health's Refseq genome. The genome data is used to create a decoy list (an important function of Salmon), allowing for greater transcriptome accuracy during the pseudoalignment. 

Transcripts will be catalogued using the GenomicFeatures and AnnotationDbi Bioconductor R packages, finding common reference names for the transcripts from the reference .gtf file. Unlike human transcript databases (which use Ensembl ID's), the yeast databases for these packages use ORF IDs, necessitating some code adjustment if using workflows relying on a symbol lookup.

These converted transcript counts are placed into a lookup dataframe with the tximport and DESeq Bioconductor packages. Log2fold change shrinkage is applied to get a clearer comparison of count changes between stages. ggplot was used to construct MA plots to get an overall picture of expression alteration distribution. pheatmap was used to view expression differences in the top 20 genes (for both early-vs-thin biofilm and early-vs-mature biofilm) by p-value across all samples. To finalize overview of the data, a Principal Component Analysis (PCA) plot was produced to compare each sample's data distribution, to see how much of the variation between samples could be attributed to larger predictor sets showing similar change patterns across the distribution.

For functional annotation, both a Gene Ontology (GO) enrichment ORA and KEGG enrichment analysis (from clusterProfiler) were performed to get a clearer functional picture of the changes in transcription activity from both a gene ontology and pathway focused analysis (respectively.)

Specific version information on software packages used may be viewed in the project readme. For the full code used to produce results, see "[Assignment_2_pseudoalignment_code.sh](https://github.com/WFeinman/BINF_6110_Assignment_2_Yeast_biofilm_transcript_analysis/blob/main/Assignment_2_pseudoalignment_code.sh "Assignment_2_pseudoalignment_code.sh")" for file download and pseudoalignment shell code, and "[Assignment_2_Functional_Annotation_Comparison.R](https://github.com/WFeinman/BINF_6110_Assignment_2_Yeast_biofilm_transcript_analysis/blob/main/Assignment_2_Functional_Annotation_Comparison.R)" for plotting and functional analysis R code.


## Results:

FastQC showed high quality for all files, with no low-quality reads recorded. This is likely due to a combination of short-read Illumina sequencing being used for the initial study, and pre-trimming of low-quality reads by the original research group already having been applied.

#### MA Plot: Thin Biofilm vs Early Biofilm
<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/0f8fe144-64c5-41b1-8832-fd20c9bd2cdd" />

#### MA Plot: Mature Biofilm vs Early Biofilm
<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/5444a583-0297-4350-a468-9741c74e6707" />

**Figures 1 and 2:** MA plots summarizing distribution of log fold change after shrinkage for transcript comparisons between biofilm development stages. Positive log fold change values indicate relative transcript increase, while negative log fold change indicates transcript decrease. Blue coloration indicates the change is statistically significant (padj < 0.05).

While both MA plots showed a similar distribution shape, the mature biofilm comparison showed a slightly higher proportion of larger log fold changes, whether for transcript increase or decrease.


<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/00763e0b-3991-4c36-b193-dcf242e421fd" />

<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/ef32890e-5be5-4634-b43d-993736bd9c1f" />

**Figures 3 and 4:** Volcano plots of transcript log fold 2 change, emphasizing display of significant transcript change rather than overall data distribution by count. These show similar results to the MA plots, though the expression level differences between stages are more easily visualized.

#### Heatmap: Top 20 Transcripts - Thin Biofilm vs Early Biofilm
<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/5d3fd7e7-fe7d-4c38-bae0-5343d519e077" />

#### Heatmap: Top 20 Transcripts - Mature Biofilm vs Early Biofilm
<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/3af6831e-c655-462d-bbc1-5d5df7309fbe" />

**Figures 5 and 6:** Heatmaps of the top 20 genes by p-value confidence displaying transcriptional differences between early biofilm and thin biofilm (Figure 5) or mature biofilm (Figure 6). Gene names are displayed by ORF-name. Cell color beside the gene names indicates increased expression (red) or decreased expression (blue). While some transcripts are shared between both heatmaps, there is a notable difference in transcript activity and the significance-listing thereof in between stages.


<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/2427c9ba-7356-4f36-8c9f-21b60b055ad9" />

**Figure 7:** PCA analysis of variation between samples. It would appear that a large number of the predictors can be grouped together for change patterns over biofilm development, as 93% of the total variance can be explained by the first 2 principle components. PC1 appears able to clearly delineate biofilm stages, while PC2 seems primarily active in distinguishing the thin biofilm stage. 


#### Enrichment Map: Transcript Activity by Linked GO terms
<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/c80b874a-888a-4352-a958-d3710e51230c" />

**Figure 8:** Enrichment map by GO-term linkage. Shows functional relations of altered gene expression in a wider biological context, grouping related terms together.


<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/1e9d955d-5e65-41be-94ba-cdf1ad143b8b" />

<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/d28f56ef-290e-4212-8de5-06489fcbfb89" />

**Figures 9 and 10:** Indicates upregulated transcript activity by gene type (Figure 9)  or metabolic process (Figure 10) in mature biofilm. Translation/Ribosome activity is distinctly elevated above other categories.


<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/efd604a6-53a5-437b-940d-31f3940a8d36" />

<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/5141eeb3-972d-4024-b319-65445f6d885e" />

**Figure 11 and 12:** Indicates downregulated transcript activity by gene type (Figure 11)  or metabolic process (Figure 12) in mature biofilm. Transmembrane transport and secondary-metabolite biosynthesis transcript activity is distinctly surpressed compared to other categories.


<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/6bd56485-819d-4190-a2ac-89d284e79852" />

<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/cd956115-0445-4e64-b67d-aa16d1404566" />

**Figures 13 and 14:** Indicates differing transcriptional activity in thin biofilm samples. Note the lower GeneRatio scale in Figure 14 compared to Figure 11.


<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/82f4ce7c-9ba9-41e9-bb11-727b0bf23909" />

<img width="1724" height="1190" alt="image" src="https://github.com/user-attachments/assets/c35ef5f1-e382-417c-81f4-c477a8aa0d91" />

**Figures 15 and 16:** Indicates slightly different metabolic priorities in thin biofilm compared to mature biofilm. While upregulated transcripts remain quite similar (Figure 15), there are a number of differences in downregulated metabolic pathway types after the top 5 most impactful, notably starch and sucrose metabolism.

## Discussion:

A large overall trend is illustrated for biofilm formation, where genes related to translation, mitochondrion organization, and translation activity are upregulated. By contrast, genes related to transmembrane transport, secondary metabolite biosynthesis, and carbon metabolism are downregulated. These are reflective of biological conditions in mature yeast flor biofilm: An emphasis on oxidative metabolism, an absence of sugar food sources, and a subsequent metabolic re-prioritization necessitating new protein architecture. (Esteve-Zarzoso et al., 2001; Fanning & Mitchell, 2012)

The most heavily upregulated significant transcript in the mature biofilm heatmap was YKR075C, a partial mRNA fragment for an uncharacterized protein. However, this protein is structurally similar to Reg1p; a binding protein within the Glc7p-Reg1p protein-phosphatase complex crucial for regulating glucose metabolism pathways. Upregulation of a glucose-metabolism inhibitor would be consistent with both the lowered observable glucose availability and the above reported downregulation of carbon metabolism pathways. As such, it is plausible that YKR075C codes for a Reg1p analogue responsible for inhibiting glucose metabolism. (Alms, 1999; Engel et al., 2022; NIH 2022)

In both heat map comparison and PCA, a distinct cluster of genes were identified specific to thin biofilm formation. One of these, YCR105W, is the most heavily upregulated gene in thin biofilm on the heatmap, but is notably downregulated in mature biofilm. YCR105W codes for ADH7; and NADP dependent alcohol dehydrogenase, involved in alcohol metabolism. A possible explanation is that ADH7 serves as a metabolic stopgap mechanism during the transition to mature biofilm, a mechanism less required when the flor yeast biofilm has fully matured.(NIH, 2026; Yadav et al., 2020)

In summary, these findings outline a number of significant transcripts active in biofilm formation, and provide a basis for future studies on wine aging or other applications of yeast metabolism.

## References:

Alms, G. R. (1999). Reg1p targets protein phosphatase 1 to dephosphorylate hexokinase II in Saccharomyces cerevisiae: Characterizing the effects of a phosphatase subunit on the yeast proteome. _The EMBO Journal_, _18_(15), 4157–4168. [https://doi.org/10.1093/emboj/18.15.4157](https://doi.org/10.1093/emboj/18.15.4157)

Esteve-Zarzoso, B., Peris-Torán, M. J., Garcı́a-Maiquez, E., Uruburu, F., & Querol, A. (2001). Yeast Population Dynamics during the Fermentation and Biological Aging of Sherry Wines. _Applied and Environmental Microbiology_, _67_(5), 2056–2061. [https://doi.org/10.1128/AEM.67.5.2056-2061.2001](https://doi.org/10.1128/AEM.67.5.2056-2061.2001)

Fanning, S., & Mitchell, A. P. (2012). Fungal Biofilms. _PLoS Pathogens_, _8_(4), e1002585. [https://doi.org/10.1371/journal.ppat.1002585](https://doi.org/10.1371/journal.ppat.1002585)

Mardanov, A. V., Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., & Ravin, N. V. (2020). Transcriptome Profile of Yeast Strain Used for Biological Wine Aging Revealed Dynamic Changes of Gene Expression in Course of Flor Development. _Frontiers in Microbiology_, _11_, 538. [https://doi.org/10.3389/fmicb.2020.00538](https://doi.org/10.3389/fmicb.2020.00538)

NCBI SRA. (2026). National Center for Biotechnology Information. [https://github.com/ncbi/sra-tools](https://github.com/ncbi/sra-tools)

NIH. (2026). ADH7 - NADP-dependent alcohol dehydrogenase. https://www.ncbi.nlm.nih.gov/datasets/gene/850469/

NIH. (2022). Saccharomyces cerevisiae S288C uncharacterized protein (YKR075C), partial mRNA. https://www.ncbi.nlm.nih.gov/nuccore/NM_001179865.1

Yadav, S., Mody, T. A., Sharma, A., & Bachhawat, A. K. (2020). A Genetic Screen To Identify Genes Influencing the Secondary Redox Couple NADPH/NADP+ in the Yeast _Saccharomyces cerevisiae_. _G3 Genes|Genomes|Genetics_, _10_(1), 371–378. [https://doi.org/10.1534/g3.119.400606](https://doi.org/10.1534/g3.119.400606)
