import os
from glob import glob
from types import SimpleNamespace
configfile: 'conf/seq-similarity-search.yaml'
config = SimpleNamespace(**config)

samples = [(lambda x: os.path.splitext(os.path.basename(x))[0])(g)
           for g in glob(f'{config.fastadir}/*.fasta')] 
index = f'{config.outdir}/{config.index_name}'
idmap = f'{index}.idmap.gz' # may not need at this time

## TODO list
# TODO add config for window size of queries
# TODO add config for shift amount of windows

rule All:
    input:
        expand(f'{config.outdir}/bed_scores/{{sample}}_scored.bed',
               sample=samples),
        expand(f'{config.outdir}/virus_beds/{{sample}}.bed', sample=samples)


rule GetModelWeights:
    output:
        'data/final.pt'
    shell:
        """
        wget -O {output} https://github.com/mchowdh200/AsMac/blob/main/model/final.pt?raw=true
        """

rule MakeFastaWindows:
    """
    Take a fasta, and create a tab sep file with format
    start  end  sequence.
    """
    input:
        f'{config.fastadir}/{{sample}}.fasta'
    output:
        f'{config.outdir}/virus_beds/{{sample}}.bed'
    shell:
        'python scripts/make_fasta_windows.py --fasta {input} --window 150 --step 50 > {output}'

# TODO
# make rule to just compute the embedding vectors of the virus sequences
# and store them in a bed file with region -> embedding vector
# TODO
# then redo the QueryIndex rule
# TODO 
# Also do a rule to compute true positives and false positives
# * query the index with a virus sequence kmer (ie embedding)
# * TP: query contains the BWA top hit
#   NOTE: bwa only give top hit, so n-nearest neighbors will necessearily
#         return results that are FP
# * FP: given the above, we will define a FP as a query did not return the BWA top hit
# * FN: with the above formulation, FP=FN
# * we can sweep the value of n-nearest neighbors to get a ROC curve





rule QueryIndex:
    """
    Using windowed sequence (bed) query the index and
    get the K-nearest neighbors.  For each window in
    the bed, compute the average(TODO is this the best way?)
    similarity of those nearest neighbors and add as column in output bed.

    Output format: {bed region}\t{average score}
    """
    input:
        index = index,
        beds = expand(f'{config.outdir}/virus_beds/{{sample}}.bed',
                      sample=samples)
    output:
        expand(f'{config.outdir}/bed_scores/{{sample}}_scored.bed',
               sample=samples)
    conda:
        'envs/asmac.yaml'
    threads:
        workflow.cores
    shell:
        # need to do some testing with k
        f"""
        python scripts/query_index.py \\
        --model-weights data/final.pt \\
        --index {{input.index}} \\
        --batch-size 64 \\
        --num-processes {{threads}} \\
        --beds {{input.beds}} \\
        --outdir {config.outdir}/bed_scores \\
        -k 10 \\
        """
        
