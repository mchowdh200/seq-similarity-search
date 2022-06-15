import sys
import numpy as np
import pandas as pd

left = pd.read_csv(sys.argv[1], sep='\t', names=['contig', 'start', 'end', 'score'])
right = pd.read_csv(sys.argv[2], sep='\t', names=['contig', 'start', 'end', 'score'])

merged = pd.merge(on=['contig', 'start', 'end'],
                  left=left, right=right, how='outer')
print(merged.head())

merged['mean_score'] = merged.apply(
    lambda x: np.nanmean([x.score_x, x.score_y]), axis=1)

merged[['contig', 'start', 'end', 'mean_score']].to_csv(
    sys.stdout, sep='\t', header=False, index=False)
