from os.path import basename, dirname, splitext
import sys
import argparse
from functools import partial
from multiprocessing import Pool

import numpy as np
import faiss
from AsMac_utility import load_pretrained, SeqIteratorDataset, makeDataLoader


parser = argparse.ArgumentParser()
parser.add_argument('--model-weights', type=str, required=True, dest='weights',
                    help='path to asmac model weights.')
parser.add_argument('--index', type=str, required=True, dest='index',
                    help='path to index.')
parser.add_argument('--batch-size', type=int, default=64, dest='batch_size',
                    help='batch size for model to operate on.')
parser.add_argument('--num-processes', type=int, default=64, dest='processes',
                    help='number of processes to use.')
parser.add_argument('--beds', type=str, required=True, nargs='+', dest='beds',
                    help='list of bed files of windowed seqs.')
parser.add_argument('-k', type=int, required=True, dest='k',
                    help='bed file of windowed seqs.')
args = parser.parse_args()


model = load_pretrained(args.weights)
index: faiss.IndexFlatL2 = faiss.read_index(args.index)
p = Pool(args.processes)

for bed in args.beds:
    dir = dirname(bed)
    sample = splitext(basename(bed))[0]
    with open(f'{dir}/{sample}_scored.bed') as out:
        dataset = SeqIteratorDataset(paths=[bed], format='bed',
                                    gzipped=False, alphabet='ATCG')
        dataloader = makeDataLoader(dataset, batch_size=args.batch_size)
        for batch in dataloader:
            embeddings = np.array(p.map(partial(model.test_embed, asnumpy=True),
                            [b['seq'] for b in batch]))
            D, _, = index.search(embeddings, k=args.k)
            for i, b in enumerate(batch):
                out.write(f"{b['interval']}\t{np.mean(D[i])}\n")
