#!/bin/bash
# Main snakemake job:
# Submits snakemake jobs as individual slurm jobs
# =============================================================================
#SBATCH -p short
#SBATCH --job-name=query-seqs
#SBATCH --mail-type=None # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=much8161@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=1gb
#SBATCH --time=1:00:00 # Time limit hrs:min:sec
#SBATCH --output=/Users/much8161/Repositories/seq-similarity-search/log/query-seqs.out
#SBATCH --error=/Users/much8161/Repositories/seq-similarity-search/log/query-seqs.err
# =============================================================================
snakemake -s query-sequences.smk \
          --use-conda --conda-frontend mamba \
          --cluster-config conf/cluster_config.yaml \
          --latency-wait 5 \
          --cluster "sbatch -J {cluster.job-name} \\
                            -t {cluster.time} \\
                            -N {cluster.nodes} \\
                            -n {cluster.ntasks} \\
                            -p {cluster.partition} \\
                            --mem={cluster.mem} \\
                            --output {cluster.output} \\
                            --error {cluster.error}" \
          --jobs 4 --cores 64