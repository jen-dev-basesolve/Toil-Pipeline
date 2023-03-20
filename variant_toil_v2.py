import os
import tempfile

from toil.common import Toil
from toil.job import Job

# Downloading Reference Genome of ecoli
def ref_genome():
    comms = "mkdir dc_workshop_1 ; cd dc_workshop_1 ; mkdir -p data/ref_genome "        # Working
    os.system(comms)
    comms = "cd dc_workshop_1 ; curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz"
    os.system(comms)
    comms = "cd dc_workshop_1 ; gunzip data/ref_genome/ecoli_rel606.fasta.gz"
    os.system(comms)


# Downloading trimmed FastQ files for faster operations
def fastq_files():
    comms = "cd dc_workshop_1 ; curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248 ; tar xvf sub.tar.gz ; mkdir data/trimmed_fastq_small ; mv sub/ data/trimmed_fastq_small"
    os.system(comms)
    comms = "cd dc_workshop_1 ; mkdir -p results/sam results/bam results/bcf results/vcf"
    os.system(comms)


# Indexing reference fasta file
def index_ref():
    print('hi')
    comms = "singularity exec bwa_latest.sif bwa index dc_workshop_1/data/ref_genome/ecoli_rel606.fasta"
    os.system(comms)


# Aligning sample sequences to reference genome
def align_fastq():
    comms = ("singularity exec bwa_latest.sif bwa mem dc_workshop_1/data/ref_genome/ecoli_rel606.fasta dc_workshop_1/data/trimmed_fastq_small/sub/SRR2584866_1.trim.sub.fastq dc_workshop_1/data/trimmed_fastq_small/sub/SRR2584866_2.trim.sub.fastq > dc_workshop_1/results/sam/SRR2584866.aligned.sam")
    os.system(comms)


# Converting SAM file to BAM file using view option in samtools
# Sorting bam files
def convrt_sort():
    comms = "singularity exec samtools_latest.sif samtools view -S -b dc_workshop_1/results/sam/SRR2584866.aligned.sam > dc_workshop_1/results/bam/SRR2584866.aligned.bam"
    os.system(comms)
    comms = "singularity exec samtools_latest.sif samtools sort -o dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam dc_workshop_1/results/bam/SRR2584866.aligned.bam"
    os.system(comms)


# Variant calling
def variant_call():
    comms = "singularity exec bcftools_latest.sif bcftools mpileup -O b -o dc_workshop_1/results/bcf/SRR2584866_raw.bcf -f dc_workshop_1/data/ref_genome/ecoli_rel606.fasta dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam"
    os.system(comms) 
    # Detecting Single Nucleotide Variants(SNV) from VCF file.
    # ploidy - number of chromosome sets in nucleus
    comms = "singularity exec bcftools_latest.sif bcftools call --ploidy 1 -m -v -o dc_workshop_1/results/vcf/SRR2584866_variants.vcf dc_workshop_1/results/bcf/SRR2584866_raw.bcf"
    os.system(comms)
    # Filter and report the SNV variants in variant calling format (VCF)
    comms = "singularity exec bcftools_latest.sif vcfutils.pl varFilter dc_workshop_1/results/vcf/SRR2584866_variants.vcf  > dc_workshop_1/results/vcf/SRR2584866_final_variants.vcf"
    os.system(comms)


# Indexing for visualization
def index_fastq():
    comms = "singularity exec samtools_latest.sif samtools index dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam"
    os.system(comms)


# Visualising SNV using t_view
def t_view():
  # Viewing with tview
  comm = "singularity exec samtools_latest.sif samtools tview dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam dc_workshop_1/data/ref_genome/ecoli_rel606.fasta"
  os.system(comm)

if __name__=='__main__':
    jobstore: str = tempfile.mkdtemp("dna_proc")
    os.rmdir(jobstore)
    options = Job.Runner.getDefaultOptions(jobstore)
    
    refgenome_job = Job(ref_genome())      # when we call the function over here, I guess it runs directly without involving toil. And when function is not called, no files are created. Need to find a way around.
    fastqfile_job = Job(fastq_files())
    indexref_job = Job(index_ref())
    refgenome_job.addChild(indexref_job)
    alignfastq_job = Job(align_fastq())
    fastqfile_job.addChild(alignfastq_job)
    convrtsort_job = Job(convrt_sort())
    alignfastq_job.addChild(convrtsort_job)
    variantcall_job = Job(variant_call())
    convrtsort_job.addChild(variantcall_job)
    indexfastq_job = Job(index_fastq())
    variantcall_job.addChild(indexfastq_job)
    tview_job = Job(t_view())
    indexfastq_job.addChild(tview_job)

    with Toil(options) as toil:
        toil.start(refgenome_job)
        toil.start(fastqfile_job)



