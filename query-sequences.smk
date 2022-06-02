import os
from glob import glob
from types import SimpleNamespace
configfile: 'conf/seq-similarity-search.yaml'
config = SimpleNamespace(**config)

samples = [(lambda x: os.path.splitext(os.bath.basename(x))[0])(g)
           for g in glob(config.fastadir)] 
index = f'{config.outdir}/{config.index_name}'
idmap = f'{index}.idmap.gz' # may not need at this time

## TODO list
# TODO add config for window size of queries
# TODO add config for shift amount of windows
# TODO add different format for iterator in AsMac for simple seq list


rule All:
    input:
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
        f'{config.fastadir}/{{virus}}.fasta'
    output:
        f'{config.outdir}/virus_beds/{{virus}}.bed'
    conda:
        'envs/asmac.yaml'
    shell:
        'python scripts/make_fasta_windows.py --fasta {input} --window 150 --step 50 > {output}'
