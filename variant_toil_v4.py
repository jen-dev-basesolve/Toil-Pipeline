from toil.common import Toil
from toil.job import Job
import os
import tempfile 
import requests
import shlex
import subprocess

def main(job,memory="2G", cores=2, disk="3G"):
    ref_job = job.addChildJobFn(ref_genome)
    indexRef_job = ref_job.addChildJobFn(index_ref,ref_job.rv())
    fastq_job = ref_job.addFollowOnJobFn(fastq_files,ref_job.rv())


def ref_genome(job):
    parent_dir = '/home/bioinfo/test_toil/variant_analysis'
    child_dir = 'dc_workshop_1/data/ref_genome'
    path = os.path.join(parent_dir,child_dir)
    os.makedirs(path,exist_ok=True)
    # Downloading fasta file
    url = "http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz"
    local_filename = 'ecoli_rel606.fasta.gz'
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(f"{path}/{local_filename}", 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                f.write(chunk)
    # extracting fasta file
    comms = ("gunzip %s/%s"%(path,local_filename))
    os.popen(comms).read()
    return parent_dir

# Indexing reference fasta file
def index_ref(job,path,memory="2G", cores=2, disk="3G"):
    comms = "singularity exec %s/bwa_latest.sif bwa index %s/dc_workshop_1/data/ref_genome/ecoli_rel606.fasta"%(path,path)
    process = os.popen(comms)
    process.read()
    job.log("".join(["Following is the path >>>>",str(path)]))


# Downloading trimmed FastQ files for faster operations
def fastq_files(job,parent_dir,memory="4G", cores=4, disk="5G"):
    child_dir = 'dc_workshop_1/data/trimmed_fastq_small'            
    path = os.path.join(parent_dir,child_dir)
    os.makedirs(path,exist_ok=True)
    
    comms = ("curl -L -o %s/sub.tar.gz https://ndownloader.figshare.com/files/14418248"%path)       # tested, runs
    os.popen(comms).read()

    comms = ("tar xvf %s/sub.tar.gz -C %s"%(path,path))     # tested, runs
    args = shlex.split(comms)
    subprocess.run(args,shell=False,capture_output=True)

    comms = "cd %s/dc_workshop_1 ; mkdir -p results/sam results/bam results/bcf results/vcf"%parent_dir     # to be tested in this function
    os.popen(comms).read()     
    
    job.addChildJobFn(align_fastq,parent_dir)


# Aligning sample sequences to reference genome
def align_fastq(job,path,memory="2G", cores=2, disk="3G"):
    comms = "cd %s ; singularity exec bwa_latest.sif bwa mem dc_workshop_1/data/ref_genome/ecoli_rel606.fasta dc_workshop_1/data/trimmed_fastq_small/sub/SRR2584866_1.trim.sub.fastq dc_workshop_1/data/trimmed_fastq_small/sub/SRR2584866_2.trim.sub.fastq > dc_workshop_1/results/sam/SRR2584866.aligned.sam"%path
    process = os.popen(comms)
    process.read()
    job.addChildJobFn(convrt_sort,path)

# Converting SAM file to BAM file using view option in samtools
# Sorting bam files
def convrt_sort(job,path,memory="2G", cores=2, disk="3G"):
    # convert sam file to bam
    comms = f"singularity exec {path}/samtools_latest.sif samtools view -S -b {path}/dc_workshop_1/results/sam/SRR2584866.aligned.sam > {path}/dc_workshop_1/results/bam/SRR2584866.aligned.bam"
    print(comms)
    os.popen(comms).read()
    # sort bam file
    comms = f"singularity exec {path}/samtools_latest.sif samtools sort -o {path}/dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam {path}/dc_workshop_1/results/bam/SRR2584866.aligned.bam"
    print(comms)
    os.popen(comms).read()
    # runs, tested.
    job.addChildJobFn(variant_call,path)

# Variant calling
def variant_call(job,path,memory="2G", cores=2, disk="3G"):
    comms = f"singularity exec {path}/bcftools_latest.sif bcftools mpileup -O b -o {path}/dc_workshop_1/results/bcf/SRR2584866_raw.bcf -f {path}/dc_workshop_1/data/ref_genome/ecoli_rel606.fasta {path}/dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam"
    print(comms)
    os.popen(comms).read()      # works, tested.
    
    # Detecting Single Nucleotide Variants(SNV) from VCF file.
    # ploidy - number of chromosome sets in nucleus
    comms = f"singularity exec {path}/bcftools_latest.sif bcftools call --ploidy 1 -m -v -o {path}/dc_workshop_1/results/vcf/SRR2584866_variants.vcf {path}/dc_workshop_1/results/bcf/SRR2584866_raw.bcf"
    os.popen(comms).read()
    # Filter and report the SNV variants in variant calling format (VCF)
    comms = f"singularity exec {path}/bcftools_latest.sif vcfutils.pl varFilter {path}/dc_workshop_1/results/vcf/SRR2584866_variants.vcf  > {path}/dc_workshop_1/results/vcf/SRR2584866_final_variants.vcf"
    os.popen(comms).read()
    job.addChildJobFn(index_fastq,path)


# Indexing for visualization
def index_fastq(job,path,memory="2G", cores=2, disk="3G"):
    comms = f"singularity exec {path}/samtools_latest.sif samtools index {path}/dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam"
    os.popen(comms).read()
    

# Visualising SNV using t_view
# def t_view():
#   path = '/home/bioinfo/singularity/variant_analysis'
#   comm = f"singularity exec {path}/samtools_latest.sif samtools tview {path}/dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam {path}/dc_workshop_1/data/ref_genome/ecoli_rel606.fasta"
#   os.popen(comm).read()


if __name__=='__main__':
    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel='INFO'
    options.clean='always'

    j1 = Job.wrapJobFn(main)

    with Toil(options) as toil:
        if  not toil.options.restart:
            print(toil.start(j1))
        else:
            print(toil.restart())
        