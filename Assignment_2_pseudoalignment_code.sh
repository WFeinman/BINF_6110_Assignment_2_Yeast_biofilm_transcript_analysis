#!/bin/bash

#Fastq tile downloads:
fasterq-dump SRR10551665
fasterq-dump SRR10551663
fasterq-dump SRR10551663

fasterq-dump SRR10551662
fasterq-dump SRR10551661
fasterq-dump SRR10551660

fasterq-dump SRR10551659
fasterq-dump SRR10551658
fasterq-dump SRR10551657

#Reference file downloads:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz

#Build Salmon index from reference genome and reference transcriptome. -p specifies threads
grep "^>" <(cat GCF_000146045.2_R64_genomic.fna) | cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

cat  GCF_000146045.2_R64_rna.fna GCF_000146045.2_R64_genomic.fna > gentrome.fa.gz

salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode


#salmon mapping for each fastq file, produces a quant.sf file in subfolder with read counts in tab seperated form
for SRR in SRR*.fastq; do
	salmon quant -i salmon_index \
	-l A -r $(echo $SRR) --validateMappings \
	-o $(echo $SRR)_transcripts_quant
done

mkdir quant_output

for SRR in SRR*_transcripts_quant; do
	
	mv $(echo $SRR) quant_output/$(echo "$SRR" |trimmed_name= sed 's/.fastq_transcripts_quant//')
done


