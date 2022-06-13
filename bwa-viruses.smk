from glob import glob
from os.path import basename, splitext
from types import SimpleNamespace

## Setup
# =============================================================================

configfile: 'conf/seq-similarity-search.yaml'
config = SimpleNamespace(**config)

# the virus strains (originally stored in fasta)
samples = [(lambda x: splitext(basename(x))[0])(g)
           for g in glob(f'{config.fastadir}/*.fasta')] 

# fastq sequences will be used as index for bwa (converted to fasta)
seqs = [basename(g).split('.', 1)[0]
        for g in glob(f'{config.fastqdir}/*.fq.gz')]

## Functions
# =============================================================================

## Rules
# =============================================================================
rule All:
    ## TODO
    input:
        expand(f'{config.outdir}/bwa_queries/{{sample}}-{{seq}}.bam',
               seq=seqs, sample=samples)

rule MakeFastaWindows:
    """
    Take a fasta, and create a new fasta with reads separated into windows.
    The read headers will contain the genomic positions
    """
    input:
        f'{config.fastadir}/{{sample}}.fasta'
    output:
        f'{config.outdir}/fasta_windows/{{sample}}-windowed.fasta'
    threads:
        1
    shell:
        """
        python scripts/make_fasta_windows.py \\
        --output-format fasta \\
        --fasta {input} \\
        --window 150 --step 50 > {output}
        """


rule Fastq2Fasta:
    input:
        f'{config.fastqdir}/{{seq}}.fq.gz'
    output:
        f'{config.outdir}/index_fasta/{{seq}}.fa.gz'
    threads:
        1
    shell:
        'bash scripts/fastq2fasta.sh {input} {output}'

rule BwaIndex:
    input: rules.Fastq2Fasta.output
    output:
        f'{config.outdir}/index_fasta/{{seq}}.fa.gz.amb',
        f'{config.outdir}/index_fasta/{{seq}}.fa.gz.ann',
        f'{config.outdir}/index_fasta/{{seq}}.fa.gz.bwt',
        f'{config.outdir}/index_fasta/{{seq}}.fa.gz.pac',
        f'{config.outdir}/index_fasta/{{seq}}.fa.gz.sa',
    threads:
        1
    conda:
        'envs/asmac.yaml'
    shell:
        'bwa index {input}'

rule BwaMem:
    input:
        index = rules.BwaIndex.output, # not explicitly used, but needed
        ref = rules.Fastq2Fasta.output,
        query = rules.MakeFastaWindows.output,
    output:
        f'{config.outdir}/bwa_queries/{{sample}}-{{seq}}.bam'
    threads:
        workflow.cores
    conda:
        'envs/asmac.yaml'
    shell:
        'bwa mem {input.ref} {input.query} -t {threads} | samtools view -bS - > {output}'
        
        
