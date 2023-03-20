import subprocess
import os
import click

exit_code = subprocess.call('./variant_steps.sh')
print(exit_code)

def statistics():
  # Getting statistics for the sorted bam file using samtools
  comm = "singularity exec samtools_latest.sif samtools flagstat dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam"
  os.system(comm)

def variants():
  comm = "less -S dc_workshop_1/results/vcf/SRR2584866_final_variants.vcf"
  os.system(comm)

def t_view():
  # Viewing with tview
  comm = "singularity exec samtools_latest.sif samtools tview dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam dc_workshop_1/data/ref_genome/ecoli_rel606.fasta"
  os.system(comm)


@click.command()
@click.option('-s', is_flag=True ,help='View statistics of the bam file')
@click.option('-v', is_flag=True ,help='View variants in VCF file')
@click.option('-t', is_flag=True ,help='T-view for viewing alignments')
def main(s,v,t):
  if s:
    statistics()
  if v:
    variants()
  if t:
    t_view()


if __name__=='__main__':
  main()
