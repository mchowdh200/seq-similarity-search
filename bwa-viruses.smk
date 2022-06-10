from glob import glob
from os.path import basename, splitext
from types import SimpleNamespace

## Setup
# =============================================================================

configfile: 'conf/seq-similarity-search.yaml'
config = SimpleNamespace(**config)

samples = [(lambda x: splitext(basename(x))[0])(g)
           for g in glob(f'{config.fastadir}/*.fasta')] 

## Functions
# =============================================================================

def seqname(w):
    """
    get the name of the seq upto the first '.'
    Keep the dirname so that the path will remain the same
    """
    return w.seq.split('.')[0]

## Rules
# =============================================================================
rule All:
    input:
        # TODO


rule MakeFastaWindows:
    """
    Take a fasta, and create a tab sep file with format
    start  end  sequence.
    """
    input:
        f'{config.fastadir}/{{sample}}.fasta'
    output:
        f'{config.outdir}/virus_beds/{{sample}}.bed'
    threads:
        1
    shell:
        'python scripts/make_fasta_windows.py --fasta {input} --window 150 --step 50 > {output}'


rule Fastq2Fasta:
    input:
        '{seq}' # full path to seq
    output:
        lambda w: f'{seqname(w)}.fa.gz'
    threads:
        1
    shell:
        'bash scripts/fastq2fasta.sh {input} {output}'

rule BwaIndex:
    input: rules.Fastq2Fasta.output
    output:
        lambda w: f'{seqname(w)}.fa.gz.amb',
        lambda w: f'{seqname(w)}.fa.gz.ann',
        lambda w: f'{seqname(w)}.fa.gz.bwt',
        lambda w: f'{seqname(w)}.fa.gz.pac',
        lambda w: f'{seqname(w)}.fa.gz.sa',
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
        query = f'{conf.fastadir}/{{sample}}.fasta',
    output:
        lambda w: f'{conf.outdir}/bwa_queries/{w.sample}-{seqname(w)}.bam'
    threads:
        workflow.cores
    conda:
        'envs/asmac.yaml'
    shell:
        'bwa mem {input.ref} {input.query} -t {threads} | samtools view -b > {output}'
        
        
