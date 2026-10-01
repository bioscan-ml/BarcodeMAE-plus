"""
Build the <name>_labels.csv files that barcodebert/datasets.py expects next to
each ITS-5M fasta file (trainset.fasta, trainset_valid.fasta, test1.fasta,
test2.fasta, test3.fasta).

Neither the Zenodo archive nor the MycoAI GitHub repo ships these CSVs --
they're derived from the raw fasta headers using mycoai's own UNITE-header
parser and taxon encoder. This reproduces the (currently commented-out)
encoding path in datasets.py's DNADataset ITS-5M loader:

    fungi_data = Data(file_path, allow_duplicates=True/False)
    tax_encoder = TaxonEncoder(data=fungi_data)   # fit ONCE, on trainset.fasta only
    labels.append(tax_encoder.encode(row))        # then reused for every other file

The encoder is fit only on trainset.fasta so its per-level vocabulary is the
training vocabulary; encoding test1/test2/test3/trainset_valid with that same
encoder is what makes species/genera absent from training collapse to the
9999999 "unknown" sentinel (mycoai.utils.UNKNOWN_INT) instead of getting their
own index -- this is intentional, not a bug, and must NOT be refit per file.

Usage:
    python build_its_labels.py --data-dir ./data/ITS-5M
"""

import argparse
import os

import pandas as pd
from mycoai.data import Data
from mycoai.data.encoders import TaxonEncoder
from mycoai.utils import LEVELS

TRAIN_FILE = "trainset.fasta"
OTHER_FILES = ["trainset_valid.fasta", "test1.fasta", "test2.fasta", "test3.fasta"]


def export_labels(data_obj, out_path, tax_encoder):
    # encode_fast is the vectorized form of per-row .encode(): same lvl_encoders,
    # same isin/transform logic, same 9999999 fallback for unseen labels -- it just
    # skips updating the inference-matrix bookkeeping, which we don't export anyway.
    # At 5.2M rows, per-row .encode() in a Python loop is impractically slow.
    encoding = tax_encoder.encode_fast(data_obj).numpy()
    out_df = pd.DataFrame(encoding, columns=LEVELS)
    out_df.insert(0, "id", data_obj.data["id"].to_numpy())
    out_df.to_csv(out_path, index=False)
    print(f"  Wrote {len(out_df)} rows -> {out_path}")


def run(data_dir):
    train_path = os.path.join(data_dir, TRAIN_FILE)
    print(f"Loading training set (fits the taxon vocabulary): {train_path}")
    train_data = Data(train_path, allow_duplicates=True)

    print("Fitting TaxonEncoder on the training set only ...")
    tax_encoder = TaxonEncoder(data=train_data)

    export_labels(train_data, train_path.replace(".fasta", "_labels.csv"), tax_encoder)
    tax_encoder.finish_training()  # finalizes inference matrices; doesn't change label ints

    for fname in OTHER_FILES:
        fpath = os.path.join(data_dir, fname)
        if not os.path.isfile(fpath):
            print(f"Skipping {fname} (not found in {data_dir})")
            continue
        print(f"Loading {fname} ...")
        allow_duplicates = "train" in fname
        data_obj = Data(fpath, allow_duplicates=allow_duplicates)
        export_labels(data_obj, fpath.replace(".fasta", "_labels.csv"), tax_encoder)


def get_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--data-dir",
        "--data_dir",
        dest="data_dir",
        required=True,
        help="Path to ITS-5M data directory (containing trainset.fasta, test1-3.fasta, trainset_valid.fasta).",
    )
    return p


def cli():
    args = get_parser().parse_args()
    run(args.data_dir)


if __name__ == "__main__":
    cli()
