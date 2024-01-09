#!/usr/bin/env python3

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "alejandrogzi"
__version__ = "0.0.2"

# -----------------------------------------------------

import numpy as np
import pandas as pd
import argparse
from gc_by_pos import format_timecourse

# -----------------------------------------------------

COLS = ["Position", "Gene", "Annotation", "codon_pos", "ref", "alt"]


def main():
    args = parse_args()
    data = format_timecourse(args.timecourse)

    if args.pos:
        data = data[data["codon_pos"] == args.pos]

    if args.gene:
        print(f"Calculating stats for gene: {args.gene}")
        data = data[data["Gene"] == args.gene]
        if args.pos:
            print(f"Calculating stats for codon position: {args.pos}")
            data = data[data["codon_pos"] == args.pos]
            print(f"Number of mutations: {len(data)} (codon position {args.pos})")
        else:
            print(
                f"Number of mutations: {len(data)} for gene {args.gene} (all codon positions)"
            )
    else:
        print(f"Number of mutations: {len(data)}")
        print(f"Type of mutations:\n{data['Annotation'].value_counts()}")

    idx = [
        col.split(":")[-1].strip() if col.startswith("AC:") else col
        for col in data.columns
    ]
    data.columns = idx


    print(f"Considering generations from {args.start} to {args.end}...\n")
    idx_map = slicer(args.start, args.end, idx[3:-5])
    data = data[idx_map]

    # Calculates the number of generations in which a mutation occurs
    muts = pd.concat(
        [
            data.iloc[:, :6],
            data.iloc[:, 6:][data.iloc[:, 6:] == 0].fillna(1).sum(axis=1),
        ],
        axis=1,
    )
    muts = muts.rename(columns={0: "muts"})
    muts = muts[
        muts["muts"] > 5
    ]  # Filter out mutations that occur in less than 5 generations
    print(f"Number of mutations that occur in more than 5 generations: {len(muts)}\n")
    print(f"Codon distribution:\n{muts['codon_pos'].value_counts()}\n")

    count = {}
    for k, v in list(zip(muts["ref"], muts["alt"])):
        count[f"{k}->{v}"] = count.get(f"{k}->{v}", 0) + 1

    print(f"Mutation distribution:\n{count}")

    return None


def slicer(start: int, end: int, idx: list) -> list:
    """Return a list of indices between start and end"""
    for k, v in enumerate(idx):
        if int(v) == start:
            st = k
            continue
        if int(v) == end:
            en = k
            break
    return COLS + idx[st : k + 1]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments"""

    parser = argparse.ArgumentParser(
        description="Calculate the GC content of a codon position based on the number of mutations that occur in that position in a given timecourse."
    )

    parser.add_argument(
        "-t",
        "--timecourse",
        metavar="timecourse",
        type=str,
        help="Timecourse file",
        required=True,
    )

    parser.add_argument(
        "-p",
        "--pos",
        metavar="pos",
        type=int,
        help="Codon position to calculate the GC content for",
        required=False,
        default=False,
    )

    parser.add_argument(
        "-g",
        "--gene",
        metavar="gene",
        type=str,
        help="Gene to calculate stats for",
        required=False,
    )

    parser.add_argument(
        "-s",
        "--start",
        metavar="start",
        type=int,
        help="Start generation",
        required=False,
        default=0,
    )

    parser.add_argument(
        "-e",
        "--end",
        metavar="end",
        type=int,
        help="End generation",
        required=False,
        default=60000,
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
