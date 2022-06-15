import argparse
import sys
from os.path import basename, splitext
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import intervaltree
from scipy.signal import savgol_filter

# TODO add "image" to the conda env

parser = argparse.ArgumentParser()
parser.add_argument('--scores-dir', dest='scores_dir',
                    type=str, required=True)
parser.add_argument('--scores-type', dest='scores_type',
                    type=str, default='similarity',
                    help='similarity score | bwa score')
args = parser.parse_args()

# scores_dir = sys.argv[1]
beds = glob(f'{args.scores_dir}/*.bed')

for bed in beds:

    interval_bins = intervaltree.IntervalTree()
    for i in range(0, 30_000, 50):
        interval_bins.addi(i, i+50, [])

    bed_intervals = intervaltree.IntervalTree()
    with open(bed) as b:
        for line in b:
            A = line.rstrip().split()
            s, e = map(int, A[1:3])
            v = float(A[3])
            bed_intervals.addi(s, e, v)

    for i in bed_intervals:
        ovlps = interval_bins.overlap(i)
        for o in ovlps:
            o.data.append(i.data)

    bins = sorted([(i.begin, i.end, np.mean(i.data))
                for i in interval_bins if i.data],
                key=lambda x: x[0])

    if args.scores_type == 'similarity':
        plt.plot([x[0] for x in bins],
                savgol_filter([1-np.sqrt(x[2])/2 for x in bins], 51, 3),
                label=basename(bed).split('_')[0])
    elif args.scores_type == 'bwa':
        plt.plot([x[0] for x in bins],
                 savgol_filter([x[2] for x in bins], 51, 3),
                 label=basename(bed.split('_')[0]))
    else:
        raise ValueError('scores_type must be "similarity" or "bwa"')

plt.legend()
plt.show()


