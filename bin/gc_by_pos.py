#!/usr/bin/env python3

# This script calculates the GC content of a codon position
# based on the number of mutations that occur in that position
# in a given timecourse. The script takes as input a timecourse
# file and the codon position to calculate the GC content for.
# The output is a dataframe with the GC content for each gene
# in the timecourse.
#
# Usage:
# ./gc_by_pos.py -t timecourse.txt -p 2 [for the third codon position]
# -----------------------------------------------------

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "alejandrogzi"
__version__ = "0.0.1"

# -----------------------------------------------------

import numpy as np
import pandas as pd
import argparse
from Bio import SeqIO

# -----------------------------------------------------

FASTA = "../../../supp/REL606.fasta"
GBK = "../../../supp/REL606.gbk"

# -----------------------------------------------------


def main():
    args = parse_args()
    dfs = calculate_gc_timecourse(format_timecourse(args.timecourse), args.pos)
    gc = make_gc_timecourse(dfs)


def format_timecourse(path: str) -> pd.DataFrame:
    data = pd.read_csv(path)
    genes = get_start_end_coords()

    data.columns = [col.strip() for col in data.columns]

    data["Annotation"] = [x.strip() for x in data["Annotation"]]
    data["Gene"] = [g.strip().split(":")[-1] for g in data["Gene"]]
    data = data[
        (data["Gene"] != "intergenic")
        & (data["P-value"] < 0.05)
        & (data["Annotation"] != "sv")
        & (data["Annotation"] != "indel")
    ]

    data["start"] = [
        genes[gene.strip()][0] if gene.strip() in genes.keys() else 0
        for gene in data["Gene"]
    ]
    data["relative_pos"] = (data["Position"] - 1) - (data["start"] - 1)
    data = data[data["relative_pos"] > 0]
    data["codon_pos"] = data["relative_pos"] % 3

    data["ref"] = [x.split("->")[0].strip() for x in data["Allele"]]
    data["alt"] = [x.split("->")[1].strip() for x in data["Allele"]]

    data = data.drop(
        columns=[
            "Test statistic",
            "P-value",
            "Deletion index",
            "Fold reduction",
            "Deletion P-value",
            "Duplication index",
            "Fold increase",
            "Duplication pvalue",
            "Passed?",
            "Allele",
        ]
    )

    for col in data.columns:
        if col.startswith("DP"):
            data.drop(columns=[col], inplace=True)

    data.reset_index(drop=True, inplace=True)

    return data


def get_start_end_coords() -> dict:
    genes = {}
    with open(GBK, "r") as gbk:
        gbk_iter = iter(gbk)
        for line in gbk_iter:
            if line.strip().startswith("gene"):
                coords = line.split()[-1].split("(")[-1].split(")")[0].split("..")
                start = int(coords[0])
                end = int(coords[1])
                name = None
                try:
                    while True:
                        nx = next(gbk_iter)
                        if "/gene" in nx or "/locus_tag":
                            name = nx.split("=")[1].strip().replace('"', "")
                            break
                except StopIteration:
                    pass
                if name:
                    pass
                    genes[name] = (start, end)
    return genes


def get_seqs_from_fasta() -> dict:
    seqs = {}
    genes = get_start_end_coords()
    with open(FASTA, "r") as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            ref_seq = record.seq
            for gene, coords in genes.items():
                start = coords[0]
                end = coords[1]
                seq = ref_seq[start - 1 : end]
                seqs[gene] = seq
    return seqs


def nt_by_codon_position(seq):
    counts = [{"G": 0, "C": 0, "A": 0, "T": 0} for _ in range(3)]

    # Iterating over the sequence by codon position
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        for pos in range(
            min(3, len(codon))
        ):  # Make sure to stay within the sequence length
            nucleotide = codon[pos]
            if nucleotide in "GCATgcat":
                counts[pos][nucleotide.upper()] += 1
    return counts


def build_mut_map(df: pd.DataFrame, threshold: int) -> dict:
    gene_muts = {}
    for idx, gene in enumerate(df["Gene"]):
        gene = gene.strip()
        allele = df.iloc[idx, -1]
        if gene not in gene_muts.keys():
            if allele >= threshold:
                gene_muts[gene] = [(df.loc[idx, "ref"], df.loc[idx, "alt"])]
            else:
                gene_muts[gene] = [(df.loc[idx, "ref"], df.loc[idx, "ref"])]
        else:
            if allele >= threshold:
                gene_muts[gene].append((df.loc[idx, "ref"], df.loc[idx, "alt"]))
            else:
                gene_muts[gene].append((df.loc[idx, "ref"], df.loc[idx, "ref"]))

    return gene_muts


def update_gc_based_on_muts(df: pd.DataFrame, pos: int, threshold=1) -> pd.DataFrame:
    seqs = get_seqs_from_fasta()
    track = []
    length = []
    for gene, seq in build_mut_map(df, threshold).items():
        if gene in seqs.keys():
            nt_by_pos = nt_by_codon_position(seqs[gene])[pos]
            if len(seq) > 1:
                for mut in seq:
                    ref = mut[0]
                    alt = mut[1]
                    nt_by_pos[ref] -= 1
                    nt_by_pos[alt] += 1
            else:
                ref = seq[0][0]
                alt = seq[0][1]
                nt_by_pos[ref] -= 1
                nt_by_pos[alt] += 1

            track.append(nt_by_pos["G"] + nt_by_pos["C"])
            length.append(sum(nt_by_pos.values()))
        else:
            print(f"{gene} not found")
            track.append(0)
            length.append(0)
    name = df.columns[-1].strip() + ":GC"
    df.drop_duplicates(subset=["Gene"], keep="first", inplace=True)
    df[name] = track
    df["length"] = length
    return df


def calculate_gc_timecourse(df: pd.DataFrame, pos: int) -> pd.DataFrame:
    dfs = []
    gc = pd.DataFrame()
    for idx, col in enumerate(df.columns):
        if col.startswith("AC"):
            tmp = df.iloc[:, [0, 1, 2, -2, -1, -3, idx]]
            # tmp = tmp[(tmp.iloc[:,-1] >= 1) & (tmp.loc[:,"codon_pos"] == pos)].reset_index(drop=True)
            tmp = tmp[(tmp.loc[:, "codon_pos"] == pos)].reset_index(drop=True)
            x = update_gc_based_on_muts(tmp, pos)
            dfs.append(x)

    # print(dfs)

    # for df in dfs:
    #     name = df.columns[-2].split(":")[1]
    #     gc[name, 0] = df.iloc[:,-2].sum() / df.iloc[:,-1].sum()

    return dfs


def make_gc_timecourse(dfs: list) -> pd.DataFrame:
    gc = pd.DataFrame()
    for df in dfs:
        name = df.columns[-2].split(":")[1]
        gc[name, 0] = df.iloc[:, -2].sum() / df.iloc[:, -1].sum()
    return gc


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate GC content by codon position"
    )
    parser.add_argument(
        "-t",
        "--timecourse",
        help="Timecourse file",
        type=str,
        metavar="FILE",
        required=True,
    )
    parser.add_argument(
        "-p", "--pos", help="Position in codon", type=int, metavar="INT", required=True
    )
    return parser.parse_args()
