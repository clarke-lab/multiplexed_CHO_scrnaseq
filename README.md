# sbo_scrnaseq_analysis

This repository enables the reproduction of the analysis of the following manuscript

# Dependancies


# 1. Preparate for analysis

## Download the FASTQ files from the EBI

In this analysis the standard 10x library preparation method was used for one sample (the TM-, t=0 sample), run on a single cell of the Chromium system. A multiplexing strategy for the remaining 5 samples was used and cells were captured from two lanes of the Chromium system. Two further libraries were produced following size selection of the cell label+SBO resulting in two further libraries that enable de-multiplexing of each individual sample. 

### mRNA libraries
```
mkdir -p raw_data/multiplexed_scRNAseq_data_1
mkdir -p raw_data/multiplexed_scRNAseq_data_2
mkdir -p raw_data/standard_scRNAseq_data

cp ../sbo_chromium_analysis/fastq/outs/fastq_path/HCNVLDRXY/CHO_untreated/*R*_001.fastq.gz raw_data/standard_scRNAseq_data
cp ../sbo_chromium_analysis/fastq/outs/fastq_path/HCNVLDRXY/SBO_mRNA_library_1/*R*_001.fastq.gz raw_data/multiplexed_scRNAseq_data_1
cp ../sbo_chromium_analysis/fastq/outs/fastq_path/HCNVLDRXY/SBO_mRNA_library_2/*R*_001.fastq.gz raw_data/multiplexed_scRNAseq_data_2
```
### SBO libraries
```
mkdir -p raw_data/SBO_library_1
mkdir -p raw_data/SBO_library_2

cp /mnt/HDD2/colin/sbo_chromium_analysis/sbo_analysis/library_1_R1.fastq.gz raw_data/SBO_library_1
cp /mnt/HDD2/colin/sbo_chromium_analysis/sbo_analysis/library_1_R2.fastq.gz raw_data/SBO_library_1

cp /mnt/HDD2/colin/sbo_chromium_analysis/sbo_analysis/library_2_R1.fastq.gz raw_data/SBO_library_2
cp /mnt/HDD2/colin/sbo_chromium_analysis/sbo_analysis/library_2_R2.fastq.gz raw_data/SBO_library_2
```


## Prepare the reference genome

At present the most complete genome assembly, PICRH, does not have an associated mtDNA sequence. We utilse the mtDNA sequence from another assembly

### Download the required filtes
```
mkdir -p reference_genome/

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna.gz \
-P reference_genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf.gz \
-P reference_genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/135/GCF_000223135.1_CriGri_1.0/GCF_000223135.1_CriGri_1.0_genomic.fna.gz \
-P reference_genome/

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/135/GCF_000223135.1_CriGri_1.0/GCF_000223135.1_CriGri_1.0_genomic.gtf.gz \
-P reference_genome/



gunzip reference_genome/*
```

### Extract the mtDNA sequence and annotation from the CriGri1.0
```
samtools faidx reference_genome/GCF_000223135.1_CriGri_1.0_genomic.fna NC_007936.1 > \
reference_genome/Crigri_1.0.mtDNA.fa

cat reference_genome/GCF_000223135.1_CriGri_1.0_genomic.gtf | \
grep ^NC_007 | grep  'protein_coding' > reference_genome/Crigri_1.0.mtDNA.gtf
```

### Merge PICRH with the mtDNA sequence and annotation
```
cat reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
reference_genome/Crigri_1.0.mtDNA.fa > \
reference_genome/reference_genome.fa

cat reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf \
reference_genome/Crigri_1.0.mtDNA.gtf > \
reference_genome/reference_genome.gtf

# Delete unnecessary files
rm reference_genome/Cri*  reference_genome/GCF*

# retain only protein coding genes 
grep -E 'gene_biotype "protein_coding"' reference_genome/reference_genome.gtf > reference_genome/reference_genome.pc.gtf 


```

# Download the ensembl genomes for RNA velocity index 

```
mkdir -p reference_genome/ensembl/picr 
wget http://ftp.ensembl.org/pub/release-104/fasta/cricetulus_griseus_picr/dna/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa.gz \
-P reference_genome/ensembl/picr 

wget http://ftp.ensembl.org/pub/release-104/gtf/cricetulus_griseus_picr/Cricetulus_griseus_picr.CriGri-PICR.104.gtf.gz \
-P reference_genome/ensembl/picr 
gunzip reference_genome/ensembl/picr/*

mkdir -p reference_genome/ensembl/chok1
wget http://ftp.ensembl.org/pub/release-104/fasta/cricetulus_griseus_crigri/dna/Cricetulus_griseus_crigri.CriGri_1.0.dna.toplevel.fa.gz \
-P reference_genome/ensembl/chok1

wget http://ftp.ensembl.org/pub/release-104/gtf/cricetulus_griseus_crigri/Cricetulus_griseus_crigri.CriGri_1.0.104.gtf.gz \
-P reference_genome/ensembl/chok1
gunzip reference_genome/ensembl/chok1/*
```

```
GENOMEDIR=reference_genome/ensembl/chok1/
raw_reference="Cricetulus_griseus_crigri.CriGri_1.0.dna.toplevel.fa"
name="ensembl_chok1_genome"
sed '/^>/ s/ .*//' $GENOMEDIR/$raw_reference > $GENOMEDIR/$name.fa
mv $GENOMEDIR/Cricetulus_griseus_crigri.CriGri_1.0.104.gtf $GENOMEDIR/$name.gtf
rm $GENOMEDIR/$raw_reference


GENOMEDIR=reference_genome/ensembl/picr/
raw_reference="Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa"
name="ensembl_picr_genome"
sed '/^>/ s/ .*//' $GENOMEDIR/$raw_reference > $GENOMEDIR/$name.fa
mv $GENOMEDIR/Cricetulus_griseus_picr.CriGri-PICR.104.gtf $GENOMEDIR/$name.gtf
rm $GENOMEDIR/$raw_reference

mkdir -p reference_genome/ensembl/hybrid
cat reference_genome/ensembl/chok1/ensembl_chok1_genome.gtf | awk '$1 == "MT"' > reference_genome/ensembl/hybrid/chok1_ensembl_mtdna.gtf

samtools faidx reference_genome/ensembl/chok1/ensembl_chok1_genome.fa
samtools faidx reference_genome/ensembl/chok1/ensembl_chok1_genome.fa MT > reference_genome/ensembl/hybrid/chok1_ensembl_mtdna.fa

cat reference_genome/ensembl/picr/ensembl_picr_genome.fa reference_genome/ensembl/hybrid/chok1_ensembl_mtdna.fa  > reference_genome/ensembl/hybrid/picr_hybrid.fasta
grep 'protein_coding' reference_genome/ensembl/picr/ensembl_picr_genome.gtf > reference_genome/ensembl/hybrid/picr_pc.gtf 
grep 'protein_coding' reference_genome/ensembl/hybrid/chok1_ensembl_mtdna.gtf > reference_genome/ensembl/hybrid/mtdna_ensembl_pc.gtf
cat reference_genome/ensembl/hybrid/picr_pc.gtf reference_genome/ensembl/hybrid/mtdna_ensembl_pc.gtf > reference_genome/ensembl/hybrid/pcir_hybrid_pc.gtf


```

# 2. Kallisto bustools

### Create the KB index
```
mkdir kallisto_index && cd kallisto_index

kb ref \
../reference_genome/reference_genome.fa \
../reference_genome/reference_genome.pc.gtf \
-i kb_cgr \
-f1 ./picr \
-g ./cho_g2t \
--overwrite 
cd ..
```

```
cd kallisto_index
kb ref \
../reference_genome/ensembl/hybrid/picr_hybrid.fasta \
../reference_genome/ensembl/hybrid/pcir_hybrid_pc.gtf \
-i kb_cgr_vel \
-f1 ./picr_vel \
-g ./cho_g2t_velo \
--overwrite \
--workflow lamanno \
-f2 intron.fa \
-c1 cdnat2c \
-c2 intront2c
cd ..
```


### KB counting
```
mkdir -p kallisto_counts/standard_scRNAseq_data

kb count -t 32 -i kallisto_index/kb_cgr -g kallisto_index/cho_g2t -x 10XV3 -o kallisto_counts/standard_scRNAseq_data --overwrite  \
raw_data/standard_scRNAseq_data/CHO_untreated_S2_L001_R1_001.fastq.gz \
raw_data/standard_scRNAseq_data/CHO_untreated_S2_L001_R2_001.fastq.gz \
raw_data/standard_scRNAseq_data/CHO_untreated_S2_L002_R1_001.fastq.gz \
raw_data/standard_scRNAseq_data/CHO_untreated_S2_L002_R2_001.fastq.gz

kb count -t 32 -i kallisto_index/kb_cgr_vel --h5ad -g kallisto_index/cho_g2t_velo -x 10XV3 -o kallisto_counts/standard_scRNAseq_data_velo --overwrite  \
-c1 kallisto_index/cdnat2c -c2 kallisto_index/intront2c --workflow lamanno --filter bustools  \
raw_data/standard_scRNAseq_data/CHO_untreated_S2_L001_R1_001.fastq.gz \
raw_data/standard_scRNAseq_data/CHO_untreated_S2_L001_R2_001.fastq.gz \
raw_data/standard_scRNAseq_data/CHO_untreated_S2_L002_R1_001.fastq.gz \
raw_data/standard_scRNAseq_data/CHO_untreated_S2_L002_R2_001.fastq.gz

mkdir -p kallisto_counts/multiplexed_scRNAseq_data_1

kb count  -t 32 -i kallisto_index/kb_cgr -g kallisto_index/cho_g2t -x 10XV3 -o kallisto_counts/multiplexed_scRNAseq_data_1 --overwrite \
raw_data/multiplexed_scRNAseq_data_1/SBO_mRNA_library_1_S3_L001_R1_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_1/SBO_mRNA_library_1_S3_L001_R2_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_1/SBO_mRNA_library_1_S3_L002_R1_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_1/SBO_mRNA_library_1_S3_L002_R2_001.fastq.gz

mkdir -p kallisto_counts/multiplexed_scRNAseq_data_1_velo

kb count  -t 32 -i kallisto_index/kb_cgr_vel -g kallisto_index/cho_g2t_velo -x 10XV3 -o kallisto_counts/multiplexed_scRNAseq_data_1_velo --overwrite \
-c1 kallisto_index/cdnat2c -c2 kallisto_index/intront2c --workflow lamanno --filter bustools  \
raw_data/multiplexed_scRNAseq_data_1/SBO_mRNA_library_1_S3_L001_R1_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_1/SBO_mRNA_library_1_S3_L001_R2_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_1/SBO_mRNA_library_1_S3_L002_R1_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_1/SBO_mRNA_library_1_S3_L002_R2_001.fastq.gz

mkdir -p kallisto_counts/multiplexed_scRNAseq_data_2

kb count -t 32 -i kallisto_index/kb_cgr -g kallisto_index/cho_g2t -x 10XV3 -o kallisto_counts/multiplexed_scRNAseq_data_2 --overwrite \
raw_data/multiplexed_scRNAseq_data_2/SBO_mRNA_library_2_S4_L001_R1_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_2/SBO_mRNA_library_2_S4_L001_R2_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_2/SBO_mRNA_library_2_S4_L002_R1_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_2/SBO_mRNA_library_2_S4_L002_R2_001.fastq.gz

mkdir -p kallisto_counts/multiplexed_scRNAseq_data_2_velo

kb count -t 32 -i kallisto_index/kb_cgr_vel -g kallisto_index/cho_g2t_velo -x 10XV3 -o kallisto_counts/multiplexed_scRNAseq_data_2_velo --overwrite \
-c1 kallisto_index/cdnat2c -c2 kallisto_index/intront2c --workflow lamanno --filter bustools  \
raw_data/multiplexed_scRNAseq_data_2/SBO_mRNA_library_2_S4_L001_R1_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_2/SBO_mRNA_library_2_S4_L001_R2_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_2/SBO_mRNA_library_2_S4_L002_R1_001.fastq.gz \
raw_data/multiplexed_scRNAseq_data_2/SBO_mRNA_library_2_S4_L002_R2_001.fastq.gz


```

## CITE-Seq SBO library

### Create white lists for CITE-Seq

```
mkdir -p CITE-seq_counts/SBO_library_1 && mkdir CITE-seq_counts/SBO_library_2

Rscript ./scripts/create_citeseq_whitelists.R
```

## SBO count
```
# library 1
CITE-seq-Count \
-R1 raw_data/SBO_library_1/library_1_R1.fastq.gz \
-R2 raw_data/SBO_library_1/library_1_R2.fastq.gz \
-t data/sbo_tags.csv \
-cells 3768 \
--start-trim 22 -cbf 1 -cbl 16 -umif 17 -umil 28 --bc_collapsing_dist 1 \
--whitelist CITE-seq_counts/SBO_library_1/whitelist.tsv \
--max-error 1 \
--sliding-window \
-o CITE-seq_counts/SBO_library_1/SBO_library_1_counts

# library 2
CITE-seq-Count \
-R1 raw_data/SBO_library_2/library_2_R1.fastq.gz \
-R2 raw_data/SBO_library_2/library_2_R2.fastq.gz \
-t data/sbo_tags.csv \
-cells 4055 \
--start-trim 22 -cbf 1 -cbl 16 -umif 17 -umil 28 --bc_collapsing_dist 1 \
--whitelist CITE-seq_counts/SBO_library_2/whitelist.tsv \
--max-error 1 \
--sliding-window \
-o CITE-seq_counts/SBO_library_2/SBO_library_2_counts
```

## Annotation
```
url=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045
wget "$url"/GCF_003668045.3_CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt \
-P reference_genome



wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_feature_table.txt.gz \
-P reference_genome
```