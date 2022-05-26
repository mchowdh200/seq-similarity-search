import argparse
import gzip

from AsMac_model import AsMac
from AsMac_utility import one_hot, SeqIteratorDataset

import faiss
import torch
from torch.utils.data import DataLoader
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
parser.add_argument('--seqs', type=str, required=True, nargs='+', dest='seqs',
                    help='list of seq file paths.')
args = parser.parse_args()

## load pretrained model
# ==============================================================================
print('loading model weights...', end='')
alphabet_size = 4 # the model was trained with ATCG only
embed_dim = 300 # number of kernel sequences
kernel_size = 20 # kernel length
net = AsMac(alphabet_size, embed_dim, kernel_size)
net_state_dict = torch.load(args.weights)
net.load_state_dict(net_state_dict)
print('done.')

## get embeddings, add to faiss index
# ==============================================================================
print('creating iterator dataset...', end='')
dataset = SeqIteratorDataset(paths=args.seqs, format=args.format,
                             gzipped=args.gzipped)
print('done.')

dataloader = DataLoader(dataset=dataset, batch_size=args.batch_size)
index = faiss.IndexFlatL2(embed_dim)
with torch.no_grad(), gzip.open(f'{args.output}.id_map.gz', 'wt') as id_map:

    print('creating faiss index...', end='')
    for records in dataloader:
        # add the embeddings (rows) to the index
        seq_oh = one_hot(records['seq'])
        embeddings = net.get_embeddings(seq_oh) \
                        .detach().numpy().astype(np.float32)
        index.add(embeddings) 
        id_map.write(f'{records["index"]}\t{records["file"]}\t{records["id"]}\n')
faiss.write_index(index, args.output)
