"""
python scripts/query_index_knn.py \
    --model-weights data/final.pt \
    --index {{input.index}} \
    --batch-size 64 \
    --num-processes {{threads}}
    --beds {{input.queries}} \
    --outdir {config.outdir}/query_result_id \
    -k 1
"""
import argparse
from functools import partial
from multiprocessing import Pool
from os.path import basename, splitext

import faiss
import numpy as np
from AsMac_utility import SeqIteratorDataset, load_pretrained, makeDataLoader
from torch import batch_norm
from torch.utils.data import dataloader


def query_index_knn(args: argparse.Namespace):
    model = load_pretrained(args.model_weights)
    index: faiss.IndexFlatL2 = faiss.read_index(args.index)
    p = Pool(args.num_processes)

    for bed in args.beds:
        sample = splitext(basename(bed))[0]
        with open(f"{args.outdir}/{sample}_ids.bed", "w") as out:
            dataset = SeqIteratorDataset(
                paths=[bed],
                format="bed",
                gzipped=False,
                alphabet="ATCG",
            )
            dataloader = makeDataLoader(dataset, batch_size=args.batch_size)

            for batch in dataloader:
                embeddings = np.array(
                    p.map(
                        partial(model.test_embed, asnumpy=True),
                        [b["seq"] for b in batch],
                    )
                )
                D, I = index.search(embeddings, k=args.k)
                for i, b in enumerate(batch):
                    out.write(f"{b['interval']}\t")
                    print(*D[i], file=out, sep=",", end='\t')
                    print(*I[i], file=out, sep=",", end='\n')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model-weights",
        type=str,
        required=True,
        dest="model_weights",
        help="Path to model weights",
    )
    parser.add_argument(
        "--index", type=str, required=True, dest="index", help="Path to index"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=64,
        dest="batch_size",
        help="Batch size for inference",
    )
    parser.add_argument(
        "--num-processes",
        type=int,
        default=1,
        dest="num_processes",
        help="Number of processes for inference",
    )
    parser.add_argument(
        "--beds",
        type=str,
        required=True,
        nargs="+",
        dest="beds",
        help="List of bed files of windowed seqs.",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=1,
        dest="k",
        help="Number of nearest neighbors to return",
    )
    args = parser.parse_args()
