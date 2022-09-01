"""
python scripts/cross_reference_bwa.py \\
        --query-bed {input.query_results} \\
        --bam1 {input.bam1} \\
        --bam2 {input.bam1} > {output}
"""
import argparse

import pysam


def query_bam(bam: str, seqname: str) -> bool:
    """
    check if the seqname is present in the bam file
    """
    for read in pysam.AlignmentFile(bam, "rb"):
        if read.query_name == seqname and (not read.is_unmapped):
            return True
    return False


def query_bam_from_seqs(seqs: list[str], bam1: str, bam2: str) -> bool:
    """
    check if any seqnames in seqs is present in the bam file
    input: seqs is a list of comma delimited pairs "bam_number,seqname"
           we only check the given bam number
    """
    for s in seqs:
        bam_num, seqname = s.split(",")
        # I hate this
        if bam_num == "1":
            bam = bam1
        elif bam_num == "2":
            bam = bam2
        else:
            raise ValueError(f"Unexpected bam number: {bam_num}")
        if query_bam(bam, seqname):
            return True
    return False


def cross_ref_bwa(args: argparse.Namespace):
    """
    check if the seqnames in the bed file are present in the bam files.
    check only aligned reads in the bam

    """
    with open(args.query_bed, "r") as f:
        for line in f:
            chr, start, end, *seq = line.rstrip().split("\t")
            print(chr, start, end, sep="\t", end="\t")

            if query_bam_from_seqs(seq, args.bam1, args.bam2):
                print("1")
            else:
                print("0")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--query-bed",
        type=str,
        required=True,
        dest="query_bed",
        help="Path to bed file with query seqname results",
    )
    parser.add_argument(
        "--bam1",
        type=str,
        required=True,
        dest="bam1",
        help="Path to bam_1 file to cross reference",
    )
    parser.add_argument(
        "--bam2",
        type=str,
        required=True,
        dest="bam2",
        help="Path to bam_2 file to cross reference",
    )
    args = parser.parse_args()
    cross_ref_bwa(args)
