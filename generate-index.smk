import os
from types import SimpleNamespace
configfile: 'conf/generate-index.yaml'
config = SimpleNamespace(**config)


rule All:
    input:
        index = f'{config.outdir}/{config.index_name}',
        id_map =  f'{config.outdir}/{config.index_name}.id_map.gz'

rule GetModelWeights:
    output:
        'data/final.pt'
    shell:
        """
        wget -O {output} https://github.com/mchowdh200/AsMac/blob/0daab7ee85d919c2fe1583c16ca824512cb220e9/model/final.pt
        """

rule IndexFastqs:
    ## TODO

rule CreateIndex:
    input:
        asmac_weights = rules.GetModelWeights.output,
        seqs = config.seqs, # list of files
    output:
        index = f'{config.outdir}/{config.index_name}',
        id_map =  f'{config.outdir}/{config.index_name}.id_map.gz'
    conda:
        'envs/asmac.yaml'
    params:
        gzipped = '--gzipped' if config.gzipped else ''
    shell:
        f"""
        python scripts/create_index.py \\
        --model-weights {{input.asmac_weights}} \\
        --output {{output.index}} \\
        --format {config.seq_filetype} \\
        {{params.gzipped}} \\
        --batch-size 64 \\
        --seqs {{input.seqs}}
        """
