#! /bin/bash

# Setting Up
# Downloading Reference Genome of ecoli

mkdir dc_workshop
cd dc_workshop
mkdir -p data/ref_genome
curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip data/ref_genome/ecoli_rel606.fasta.gz



# Downloading trimmed FastQ files for faster operations
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz
mkdir data/trimmed_fastq_small
mv sub/ data/trimmed_fastq_small
mkdir -p results/sam results/bam results/bcf results/vcf

# Running commands in containers
# Indexing reference fasta file
singularity exec bwa_latest.sif bwa index dc_workshop/data/ref_genome/ecoli_rel606.fasta

# Aligning sample sequences to reference genome
singularity exec bwa_latest.sif bwa mem dc_workshop/data/ref_genome/ecoli_rel606.fasta dc_workshop/data/trimmed_fastq_small/sub/SRR2584866_1.trim.sub.fastq dc_workshop/data/trimmed_fastq_small/sub/SRR2584866_2.trim.sub.fastq > dc_workshop/results/sam/SRR2584866.aligned.sam

# Converting SAM file to BAM file using view option in samtools
singularity exec samtools_latest.sif samtools view -S -b dc_workshop/results/sam/SRR2584866.aligned.sam > dc_workshop/results/bam/SRR2584866.aligned.bam


# Sorting bam files
singularity exec samtools_latest.sif samtools sort -o dc_workshop/results/bam/SRR2584866.aligned.sorted.bam dc_workshop/results/bam/SRR2584866.aligned.bam

# Getting statistics for the sorted bam file using samtools
singularity exec samtools_latest.sif samtools flagstat dc_workshop/results/bam/SRR2584866.aligned.sorted.bam

# Variant calling
singularity exec bcftools_latest.sif bcftools mpileup -O b -o dc_workshop/results/bcf/SRR2584866_raw.bcf -f dc_workshop/data/ref_genome/ecoli_rel606.fasta dc_workshop/results/bam/SRR2584866.aligned.sorted.bam


# Detecting Single Nucleotide Variants(SNV) from VCF file.
# ploidy - number of chromosome sets in nucleus
singularity exec bcftools_latest.sif bcftools call --ploidy 1 -m -v -o dc_workshop/results/vcf/SRR2584866_variants.vcf dc_workshop/results/bcf/SRR2584866_raw.bcf


# Filter and report the SNV variants in variant calling format (VCF)
singularity exec bcftools_latest.sif vcfutils.pl varFilter dc_workshop/results/vcf/SRR2584866_variants.vcf  > dc_workshop/results/vcf/SRR2584866_final_variants.vcf

# Exploring vcf files
less -S dc_workshop/results/vcf/SRR2584866_final_variants.vcf


# Assess the alignment (visualization)
# Indexing for visualization
singularity exec samtools_latest.sif samtools index dc_workshop/results/bam/SRR2584866.aligned.sorted.bam


# Viewing with tview
singularity exec samtools_latest.sif samtools tview dc_workshop/results/bam/SRR2584866.aligned.sorted.bam dc_workshop/data/ref_genome/ecoli_rel606.fasta






























