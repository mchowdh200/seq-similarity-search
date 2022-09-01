"""
python scripts/get_seqnames_from_idmap.py \\
    --idmap {input.idmap} \\
    --bed {input.query_results} \\
    > {output}
"""
import argparse
import gzip
from os.path import basename, splitext


def load_idmap(path: str) -> dict[str, tuple[str, str]]:
    """format of idmap: id\tfastq_path\tseqname.
    need to strip the directory and extension from the fastq_path
    """
    with gzip.open(path, "rt") as f:
        idmap: dict[str, tuple[str, str]] = {}
        for line in f:
            id, fastq_path, seqname = line.rstrip().split("\t")
            # fastq_number = splitext(basename(fastq_path))[0].split("_")[-1]
            fastq_number = basename(fastq_path).split('.')[0].split("_")[-1]
            # just keep the suffix of the fastq (1 or 2)
            idmap[id] = (fastq_number, seqname)
        return idmap


def get_seqnames(idmap: dict[str, tuple[str, str]], bed: str):
    with open(bed, "r") as f:
        for line in f:
            # interval, D, I. Tab separated
            # D and I are comma separated.
            chr, start, end, _, ids = line.rstrip().split("\t")
            ids = ids.split(",")
            seqnames = [",".join(idmap[id]) for id in ids]
            print(chr, start, end, *seqnames, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--idmap",
        type=str,
        required=True,
        dest="idmap",
        help="Path to idmap file",
    )
    parser.add_argument(
        "--bed",
        type=str,
        required=True,
        dest="bed",
        help="Path to bed file",
    )
    args = parser.parse_args()
    idmap = load_idmap(args.idmap)
    get_seqnames(idmap, args.bed)
