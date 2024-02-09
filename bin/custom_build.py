#!/usr/bin/env python3

# This script is just a compilation of the functions used
# to perform the statistical anaylsis of the mutational 
# distribution based on Lenski's series.
# Note that some of the functions here are repeated.
# -----------------------------------------------------

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "alejandrogzi"
__version__ = "0.0.1"

# -----------------------------------------------------

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
from pathlib import Path
import matplotlib.pyplot as plt
from functools import reduce
from scipy.stats import ks_2samp

# -----------------------------------------------------

GENOME_SIZE = 4629812
SIZE = 25298 // 2
ITERATIONS = 100
COLS = ["Position","Gene", "Annotation", "codon_pos", "ref", "alt", "start", "relative_pos"]
MUTATORS = ["m1", "m2", "m3", "m4", "p3", "p6"]
NONMUTATORS = ["m5", "m6", "p1", "p2", "p4", "p5"]
at_to_gc = ["A->G", "A->C", "T->G", "T->C"]
gc_to_at = ["G->A", "G->T", "C->A", "C->T"]
PATH = Path("../supp/dbs/")


def main() -> (pd.DataFrame, np.array):
    dfs = []

    for file in list(PATH.glob("*.txt")):
        pop = str(file).split("/")[-1].split("_")[0]
        rs = format_timecourse(file)
        df = rs[0][["Position", "Annotation"]]

        # only consider frequencies up to 60,000 generations (no clones)
        df["count"] = rs[2][rs[2]<1].fillna(1).loc[:,:60000].sum(axis=1)
        dfs.append(df)

    counts = pd.concat(dfs, ignore_index=True)
    counts = counts.groupby("Position", as_index=False).agg({"Annotation": "first", "count": "sum"})

    print(counts.head())
    genome_vec = np.array([0]*GENOME_SIZE)
    genome_vec[counts["Position"].to_list()] = counts["count"]

    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16,4))
    
    ax1.set_title("Mutation distribution across genome")
    ax1.set_xlim(0,4629812)
    ax1.set_xlabel("Genome pos [bp]")

    ax2.set_title("Mutation distribution across genome")
    ax2.set_xlim(0,4629812)
    ax2.set_xlabel("Genome pos [bp]")

    ax3.set_title("KS-statistic vs p-value")
    ax3.set_xlabel("p-value")
    ax3.set_ylabel("KS-statistic")


    for x in range(ITERATIONS):
        vec1 = np.random.choice(counts.index, size=SIZE, replace=False)
        set1 = counts.iloc[vec1]
        
        # generate the complementary subset for set2
        all_indices = set(counts.index)
        selected_indices = set(vec1)
        complement_indices = np.array(list(all_indices - selected_indices))
        set2 = counts.iloc[complement_indices]

        set1["Position"].plot(kind="density", ax=ax1, bw_method=0.09, alpha=.2, color="blue", linestyle="dashed")
        set2["Position"].plot(kind="density", ax=ax2, bw_method=0.09, alpha=.2, color="red", linestyle="dashed")

        # kolmogorov-smirnov test
        ks_statistic, p_value = ks_2samp(set1["Position"].sort_values(), set2["Position"].sort_values())
        ax3.scatter(p_value, ks_statistic, color="grey", alpha=.5)
        ax3.axvline(x=0.05, color="red")

    plt.tight_layout(pad=2.5)
    plt.savefig(
    "./dist_ks_sims.pdf",
    format="pdf",
    bbox_inches="tight",
    dpi=300,
)  

    return (counts, genome_vec)


# 1) extract coordinates for each gene and store them in a list of tuples
def get_start_end_coords() -> dict:
    genes = {}
    with open("../supp/REL606.gbk", "r") as gbk:
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
                            name = nx.split("=")[1].strip().replace('"', '')
                            break
                except StopIteration:
                    pass
                if name:
                    pass
                    genes[name] = (start,end)
    return genes


def get_seqs_from_fasta() -> dict:
    seqs = {}
    genes = get_start_end_coords()
    with open("../supp/REL606.fasta", "r") as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            ref_seq = record.seq
            for gene, coords in genes.items():
                start = coords[0]
                end = coords[1]
                seq = ref_seq[start-1:end]
                seqs[gene] = seq
    return seqs


def nt_by_codon_position(seq):
    counts = [{'G': 0, 'C': 0, 'A':0, 'T':0} for _ in range(3)]

    # Iterating over the sequence by codon position
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        for pos in range(min(3, len(codon))):  # Make sure to stay within the sequence length
            nucleotide = codon[pos]
            if nucleotide in 'GCATgcat':
                counts[pos][nucleotide.upper()] += 1
    return counts


def format_timecourse(file: str, essential=None) -> tuple:
    data = pd.read_csv(file)

    genes = get_start_end_coords()
    
    data.columns = [col.strip() for col in data.columns]
    
    data["Annotation"] = [x.strip() for x in data["Annotation"]]
    data["Gene"] = [g.strip().split(":")[-1] for g in data["Gene"]]
    data = data[(data["Gene"] != "intergenic") & (data["P-value"] < 0.05) & (data["Annotation"] != "sv") & (data["Annotation"] != "indel") & (data["Passed?"] == " PASS") & (data["Annotation"] != "noncoding")]
    
    data["start"] = [genes[gene.strip()][0] if gene.strip() in genes.keys() else 0 for gene in data["Gene"]]
    data["relative_pos"] = (data["Position"] - 1) - (data["start"]- 1)
    data = data[data["relative_pos"] > 0]
    data["codon_pos"] = data["relative_pos"] % 3
    
    data["ref"] = [x.split("->")[0].strip() for x in data["Allele"]]
    data["alt"] = [x.split("->")[1].strip() for x in data["Allele"]]
    
    data = data.drop(columns=['Test statistic', 'P-value','Deletion index', 'Fold reduction', 'Deletion P-value','Duplication index', 'Fold increase', 'Duplication pvalue', 'Passed?', "Allele"])

    # data = pd.merge(data, ESS_DB, on="Gene", how="left")
    # if essential:
    #     data = data[data[ESSENCIALITY] == True]         
    #
    guide = data[COLS]
    depths = data[[col for col in data.columns if col.startswith("DP")]]
    depths.columns = [int(col.split(":")[-1].strip()) for col in depths.columns]
    
    alleles = data[[col for col in data.columns if col.startswith("AC")]]
    alleles.columns = [int(col.split(":")[-1].strip()) for col in alleles.columns]

    return (guide, depths, alleles)


if __name__ == "__main__":
    main()
