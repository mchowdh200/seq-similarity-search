import argparse
from typing import Iterator

class SeqRecord:
    def __init__(self, contig, start, end, seq):
        self.contig: str = contig
        self.start: int = start
        self.end: int = end
        self.seq: str = seq
    def __repr__(self):
        return '\t'.join([self.contig, str(self.start), str(self.end), self.seq])

def sliding_window(contig: str, seq: str, window: int, step: int) -> Iterator[SeqRecord]:
    N = len(seq)
    for i in range(0, N, step):
        end = min(i+window, N)
        yield SeqRecord(contig, i, end, seq[i:end])
        if end == N:
            return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str, dest='fasta', required=True)
    parser.add_argument('--window', type=int, dest='window', required=True)
    parser.add_argument('--step', type=int, dest='step', required=True)
    args = parser.parse_args()

    with open(args.fasta) as f:
        contig = f.readline().split()[0][1:]
        seq = ''.join(i.rstrip() for i in f.readlines())
        for rec in sliding_window(contig, seq, 150, 50):
            print(rec)
