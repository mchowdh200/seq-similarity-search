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
        wget -O {output} https://github.com/mchowdh200/AsMac/blob/main/model/final.pt?raw=true
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
        mkdir -p log
        python scripts/create_index.py \\
        --model-weights {{input.asmac_weights}} \\
        --output {{output.index}} \\
        --format {config.seq_filetype} \\
        {{params.gzipped}} \\
        --batch-size 64 \\
        --seqs {{input.seqs}}
        """
