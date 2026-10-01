"""
Export the actual deduplicated, genus-level ITS test sets (Yeast and
Filamentous) as standalone FASTA files, for sharing outside this repo.

These are the specimens actually used for evaluation in the paper: queries
whose exact barcode is NOT in the training set (no leakage) and whose species
is novel (not seen in training) but whose genus is known, i.e. task ==
"genus_level" in the CSVs produced by analyze_its_overlap.py --export-dir.
Counts: 526 Yeast, 3,136 Filamentous (matches the paper).

Usage:
    python export_dedup_test_sets.py \
        --data-dir ./data/ITS-5M \
        --tasks-dir ./data/ITS-5M/tasks \
        --out-dir ./data/ITS-5M/dedup_exports
"""

import argparse
import os

TEST_SETS = [
    ("test1", "yeast"),
    ("test2", "filamentous"),
]


def read_fasta(path):
    """Return list of (header, sequence) tuples. header is the full line
    after '>', including the UNITE taxonomy string -- the leading token up to
    the first '|' is the specimen id used in the tasks.csv 'id' column."""
    records = []
    header, seq_parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def write_fasta(records, path):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def load_genus_level_ids(tasks_csv):
    import pandas as pd

    df = pd.read_csv(tasks_csv)
    return set(df.loc[df["task"] == "genus_level", "id"].astype(str))


def run(data_dir, tasks_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for tag, name in TEST_SETS:
        fasta_path = os.path.join(data_dir, f"{tag}.fasta")
        tasks_path = os.path.join(tasks_dir, f"{tag}_tasks.csv")
        out_path = os.path.join(out_dir, f"{tag}_{name}_dedup_genus_level.fasta")

        print(f"Loading {fasta_path} ...")
        records = read_fasta(fasta_path)
        keep_ids = load_genus_level_ids(tasks_path)

        kept = [(header, seq) for header, seq in records if header.split("|", 1)[0] in keep_ids]
        write_fasta(kept, out_path)
        print(f"  {len(kept)} / {len(records)} specimens kept (genus_level only) -> {out_path}")
        if len(kept) != len(keep_ids):
            print(
                f"  WARNING: {len(keep_ids)} ids expected from {tasks_path} but only "
                f"{len(kept)} matched in {fasta_path} -- check tag/file alignment."
            )


def get_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", "--data_dir", dest="data_dir", required=True, help="Dir with test1.fasta/test2.fasta.")
    p.add_argument(
        "--tasks-dir",
        "--tasks_dir",
        dest="tasks_dir",
        required=True,
        help="Dir with test1_tasks.csv/test2_tasks.csv (from analyze_its_overlap.py --export-dir).",
    )
    p.add_argument("--out-dir", "--out_dir", dest="out_dir", required=True, help="Where to write the exported fasta files.")
    return p


def cli():
    args = get_parser().parse_args()
    run(args.data_dir, args.tasks_dir, args.out_dir)


if __name__ == "__main__":
    cli()
