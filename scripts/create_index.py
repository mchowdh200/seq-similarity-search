from multiprocessing import Pool
import argparse
import gzip

from AsMac_model import AsMac
from AsMac_utility import (load_pretrained, SeqIteratorDataset,
                           makeDataLoader, formatBatchMetadata)

import faiss
import torch
import numpy as np

## Get args
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument('--model-weights', type=str, required=True, dest='weights',
                    help='path to asmac model weights.')
parser.add_argument('--output', type=str, required=True, dest='output',
                    help='output path/name for faiss index.')
parser.add_argument('--format', type=str, required=True, dest='format',
                    help='seq file format -- fast{a|q}.')
parser.add_argument('--gzipped', action='store_true', default=False,
                    dest='gzipped', help='Is the file gzipped?')
parser.add_argument('--batch-size', type=int, default=64, dest='batch_size',
                    help='batch size for model to operate on.')
parser.add_argument('--num-processes', type=int, default=64, dest='processes',
                    help='number of processes to use.')
parser.add_argument('--seqs', type=str, required=True, nargs='+', dest='seqs',
                    help='list of seq file paths.')
args = parser.parse_args()


embed_dim = 300
model = load_pretrained(args.weights)
## TODO change iterator to be case insensitive
dataset = SeqIteratorDataset(paths=args.seqs, format=args.format,
                             gzipped=args.gzipped, alphabet='ATCG')
print('done.')

dataloader = makeDataLoader(dataset, batch_size=args.batch_size,
                            num_workers=args.processes)
index = faiss.IndexFlatL2(embed_dim)
p = Pool(args.processes)
with gzip.open(f'{args.output}.id_map.gz', 'wt') as id_map:

    print('creating faiss index...', end='')
    for batch in dataloader:
        # add the embeddings (rows) to the index
        embeddings = p.map(model.test_embed, [b['seq'] for b in batch])
        for s in formatBatchMetadata(batch):
            id_map.write(s)
            id_map.write('\n')
        index.add(np.array(embeddings))
faiss.write_index(index, args.output)
