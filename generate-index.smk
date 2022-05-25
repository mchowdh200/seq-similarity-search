import os
from types import SimpleNamespace
configfile: 'conf/generate-index.yaml'
config = SimpleNamespace(**config)


rule All:
    input:
        index = f'{config.outdir}/{config.index_name}',
        id_map =  f'{config.outdir}/{config.index_name}.id_map.gz'

rule CompileAsMac:
    output:
        'AsMac/_softnw.cpython-38-x86_64-linux-gnu.so'
    conda:
        'envs/asmac.yaml'
    shell:
        """
        cd AsMac
        python setup_softnw.py build_ext --inplace
        """

rule IndexFastqs:
    ## TODO

rule CreateIndex:
    input:
        asmac_compiled = rules.CompileAsMac.output,
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
        scripts/create_index.py \\
        --model-weights AsMac/model/final.pt \\
        --output {{output.index}} \\
        --format {config.seq_filetype} \\
        {{params.gzipped}} \\
        --batch-size 64 \\
        --seqs {{input.seqs}}
        """
