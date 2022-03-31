# Multiplexed CHO single cell RNA-seq

The code contained in this repositority enable the reproduction of the results of:

Tzani *et. al* 2022. Found in translation: Microproteins are a new class of host cell impurity in mAb drug products

## Dependancies

| Software | R packages      ||
| ------------- | --------------- | --------------- |

## to be removed when made public
```
# set the cell ranger path
export PATH=/mnt/HDD2/colin/bin/cellranger-6.0.2:$PATH
 
# run cell cellranger mkfastq
# converts the bcl files based on the indexes used
cellranger mkfastq \
--id fastq \
--localcores=32 \
--run /mnt/HDD2/colin/NT_003_Novarun_14JUN21/Files/ \
--csv /mnt/HDD2/colin//multiplexed_CHO_scrnaseq/data/sbo_library_sample_sheet.csv
```
 
# cat the sbo fastq files and rename
```
mkdir -p raw_data/sbo_libraries
cat fastq/outs/fastq_path/HCNVLDRXY/SBO_tag_library_1_S17_L001_R1_001.fastq.gz \
fastq/outs/fastq_path/HCNVLDRXY/SBO_tag_library_1_S17_L002_R1_001.fastq.gz > \
raw_data/sbo_libraries/sbo_library_1_R1.fastq.gz
 
cat fastq/outs/fastq_path/HCNVLDRXY/SBO_tag_library_1_S17_L001_R2_001.fastq.gz \
fastq/outs/fastq_path/HCNVLDRXY/SBO_tag_library_1_S17_L002_R2_001.fastq.gz > \
raw_data/sbo_libraries/sbo_library_1_R2.fastq.gz
 
cat fastq/outs/fastq_path/HCNVLDRXY/SBO_tag_library_2_S18_L001_R1_001.fastq.gz \
fastq/outs/fastq_path/HCNVLDRXY/SBO_tag_library_2_S18_L002_R1_001.fastq.gz > \
raw_data/sbo_libraries/sbo_library_2_R1.fastq.gz
 
cat fastq/outs/fastq_path/HCNVLDRXY/SBO_tag_library_2_S18_L001_R2_001.fastq.gz \
fastq/outs/fastq_path/HCNVLDRXY/SBO_tag_library_2_S18_L002_R2_001.fastq.gz > \
raw_data/sbo_libraries/sbo_library_2_R2.fastq.gz
 
echo $(zcat raw_data/sbo_libraries/sbo_library_1_R2.fastq.gz |wc -l)/4|bc
echo $(zcat data/fastq/SBO_library_1/library_1_R2.fastq.gz |wc -l)/4|bc
```
 
# cp and rename the cellular RNA files
```
mkdir raw_data/mrna_libraries
cp ../sbo_chromium_analysis/fastq/outs/fastq_path/HCNVLDRXY/SBO_mRNA_library_1/*R*_001.fastq.gz raw_data/mrna_libraries
cp ../sbo_chromium_analysis/fastq/outs/fastq_path/HCNVLDRXY/SBO_mRNA_library_2/*R*_001.fastq.gz raw_data/mrna_libraries
```
 
# 1. Preparate for analysis
 
## Download the FASTQ files from the EBI
 
1. In this analysis the standard 10x library preparation method was used for one sample (the TM-, t=0 sample), run on a single cell of the Chromium system. A multiplexing strategy for the remaining 5 samples was used and cells were captured from two lanes of the Chromium system. Two further libraries were produced following size selection of the cell label+SBO resulting in two further libraries that enable de-multiplexing of each individual sample.
 
### mRNA libraries
```
```
### SBO libraries
```
 
```
 
## Prepare the reference genome
 
At present the most complete genome assembly, PICRH, does not have an associated mtDNA sequence. We utilse the mtDNA sequence from another assembly
 
### Download the required files
```
mkdir -p reference_genome/
 
# 1. Chinese hamster genome
wget http://ftp.ensembl.org/pub/release-105/gtf/cricetulus_griseus_picr/Cricetulus_griseus_picr.CriGri-PICR.105.gtf.gz -P reference_genome
 
wget http://ftp.ensembl.org/pub/release-105/fasta/cricetulus_griseus_picr/dna/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa.gz -P reference_genome
 
 
# 2. CHO K1 genome
wget http://ftp.ensembl.org/pub/release-105/fasta/cricetulus_griseus_crigri/dna/Cricetulus_griseus_crigri.CriGri_1.0.dna.toplevel.fa.gz -P reference_genome
 
wget http://ftp.ensembl.org/pub/release-105/gtf/cricetulus_griseus_crigri/Cricetulus_griseus_crigri.CriGri_1.0.105.gtf.gz -P reference_genome
 
gunzip reference_genome/*
```
 
### Extract the mtDNA sequence from the CriGri1.0 and merge with the cgr reference genome
```
samtools faidx reference_genome/Cricetulus_griseus_crigri.CriGri_1.0.dna.toplevel.fa MT > \
reference_genome/chok1_mtdna.fasta
 
cat reference_genome/Cricetulus_griseus_crigri.CriGri_1.0.105.gtf | \
grep ^MT  > reference_genome/chok1_mtdna.gtf
 
cat reference_genome/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa \
reference_genome/chok1_mtdna.fasta > \
reference_genome/combined_reference_genome.fasta
 
cat reference_genome/Cricetulus_griseus_picr.CriGri-PICR.105.gtf reference_genome/chok1_mtdna.gtf > reference_genome/combined_reference_genome.gtf
 
../bin/cellranger-6.0.2/bin/cellranger mkgtf \
reference_genome/combined_reference_genome.gtf \
reference_genome/combined_reference_genome.filtered.gtf \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:lncRNA \
--attribute=gene_biotype:IG_J_gene \
--attribute=gene_biotype:IG_V_gene \
--attribute=gene_biotype:TR_C_gene \
--attribute=gene_biotype:TR_J_gene \
--attribute=gene_biotype:TR_V_gene
 
 
# clean up
rm reference_genome/Cricetulus_griseus* 
rm reference_genome/chok1*
```
 
# 2. Kallisto bustools
 
### Create the KB index
```
mkdir kallisto_index && cd kallisto_index
 
kb ref \
../reference_genome/combined_reference_genome.fasta \
../reference_genome/combined_reference_genome.filtered.gtf \
-i kb_cgr \
-f1 ./picr \
-g ./cho_g2t \
--overwrite
cd ..
 
```
 
### KB counting
```
mkdir -p kallisto_counts/mrna_library_1
 
kb count  -t 32 -i kallisto_index/kb_cgr -g kallisto_index/cho_g2t -x 10XV3 -o kallisto_counts/mrna_library_1 --overwrite \
raw_data/mrna_libraries/SBO_mRNA_library_1_S3_L001_R1_001.fastq.gz \
raw_data/mrna_libraries/SBO_mRNA_library_1_S3_L001_R2_001.fastq.gz \
raw_data/mrna_libraries/SBO_mRNA_library_1_S3_L002_R1_001.fastq.gz \
raw_data/mrna_libraries/SBO_mRNA_library_1_S3_L002_R2_001.fastq.gz
 
mkdir -p kallisto_counts/mrna_library_2
 
kb count -t 32 -i kallisto_index/kb_cgr -g kallisto_index/cho_g2t -x 10XV3 -o kallisto_counts/mrna_library_2 --overwrite \
raw_data/mrna_libraries/SBO_mRNA_library_2_S4_L001_R1_001.fastq.gz \
raw_data/mrna_libraries/SBO_mRNA_library_2_S4_L001_R2_001.fastq.gz \
raw_data/mrna_libraries/SBO_mRNA_library_2_S4_L002_R1_001.fastq.gz \
raw_data/mrna_libraries/SBO_mRNA_library_2_S4_L002_R2_001.fastq.gz
 
 
```
 
## CITE-Seq SBO library
 
### Create white lists for CITE-Seq
 
```
mkdir -p sbo_counts/SBO_library_1 && mkdir sbo_counts/SBO_library_2
 
Rscript ./r_scripts/create_citeseq_whitelists.R
```
 
## SBO count
```
# library 1
CITE-seq-Count \
-R1 raw_data/sbo_libraries/sbo_library_1_R1.fastq.gz \
-R2 raw_data/sbo_libraries/sbo_library_1_R2.fastq.gz \
-t data/sbo_tags.csv \
-cells 5000 \
--start-trim 22 -cbf 1 -cbl 16 -umif 17 -umil 28 --bc_collapsing_dist 1 \
--whitelist sbo_counts/SBO_library_1/whitelist.tsv \
--max-error 1 \
--sliding-window \
-o sbo_counts/SBO_library_1/SBO_library_1_counts_test
 
# library 2
CITE-seq-Count \
-R1 raw_data/sbo_libraries/sbo_library_2_R1.fastq.gz \
-R2 raw_data/sbo_libraries/sbo_library_2_R2.fastq.gz \
-t data/sbo_tags.csv \
-cells 5000 \
--start-trim 22 -cbf 1 -cbl 16 -umif 17 -umil 28 --bc_collapsing_dist 1 \
--whitelist sbo_counts/SBO_library_2/whitelist.tsv \
--max-error 1 \
--sliding-window \
-o sbo_counts/SBO_library_2/SBO_library_2_counts_test
```