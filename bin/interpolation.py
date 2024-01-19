#!/usr/bin/env python3


# -----------------------------------------------------

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "alejandrogzi"
__version__ = "0.0.1"

# -----------------------------------------------------

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from functools import reduce
from scipy.interpolate import interp1d
from gc_by_pos import get_start_end_coords

# -----------------------------------------------------

FASTA = "../supp/REL606.fasta"
GBK = "../supp/REL606.gbk"
COLS = [
    "Position",
    "Gene",
    "Annotation",
    "codon_pos",
    "ref",
    "alt",
    "start",
    "relative_pos",
]
MUTATORS = ["m1", "m2", "m3", "m4", "p3", "p6"]
NONMUTATORS = ["m5", "m6", "p1", "p2", "p4", "p5"]
AT_TO_GC = ["A->G", "A->C", "T->G", "T->C"]
GC_TO_AT = ["G->A", "G->T", "C->A", "C->T"]
MIN_COVERAGE = 5
MIN_DEPTH = 5
LINE_COLOR_MAP = {
    "m5": "#4A1486",
    "p2": "#807DBA",
    "p4": "#084594",
    "p1": "#4292C6",
    "p5": "#005A32",
    "m6": "#41AB5D",
    "m1": "#8C2D04",
    "m2": "#CC4C02",
    "m3": "#B10026",
    "m4": "#E31A16",
    "p3": "#FC4E2A",
    "p6": "#FD8D3C",
}
LINESTYLE = "o-"
POS_COLORS = {0: "#084594", 1: "#4292C6", 2: "#CC4C02"}
POS_DIR_COLORS = {
    "0-AT->GC": "#4A1486",
    "1-AT->GC": "#4292C6",
    "2-AT->GC": "#005A32",
    "0-GC->AT": "#B10026",
    "1-GC->AT": "#E31A16",
    "2-GC->AT": "#FD8D3C",
}

# -----------------------------------------------------


def main():
    args = parse_args()
    path = Path(args.path)
    files = list(path.glob("*.txt"))
    mutation_trajectories, avg_timecourse = runner(files)


def format_timecourse(file: str, essential=None) -> tuple:
    data = pd.read_csv(file)

    genes = get_start_end_coords()

    data.columns = [col.strip() for col in data.columns]

    data["Annotation"] = [x.strip() for x in data["Annotation"]]
    data["Gene"] = [g.strip().split(":")[-1] for g in data["Gene"]]
    data = data[
        (data["Gene"] != "intergenic")
        & (data["P-value"] < 0.05)
        & (data["Annotation"] != "sv")
        & (data["Annotation"] != "indel")
        & (data["Passed?"] == " PASS")
        & (data["Annotation"] != "noncoding")
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

    data = pd.merge(data, ESS_DB, on="Gene", how="left")
    if essential:
        data = data[data[ESSENCIALITY] == True]

    guide = data[COLS]
    depths = data[[col for col in data.columns if col.startswith("DP")]]
    depths.columns = [int(col.split(":")[-1].strip()) for col in depths.columns]

    alleles = data[[col for col in data.columns if col.startswith("AC")]]
    alleles.columns = [int(col.split(":")[-1].strip()) for col in alleles.columns]

    return (guide, depths, alleles)


def interpolate(dfs: tuple, filter=False) -> np.array:
    guide, depths, alleles = dfs
    avg_depths = depths.mean()

    times = np.array(depths.columns)
    times_no_clones = times[times < 100000]
    Ms = np.zeros_like(times_no_clones) * 1.0

    depths *= (depths >= MIN_COVERAGE) * (avg_depths >= MIN_COVERAGE)
    alleles *= (depths >= MIN_COVERAGE) * (avg_depths >= MIN_COVERAGE)  # *(alleles>=2)

    depths = depths.reset_index(drop=True)
    alleles = alleles.reset_index(drop=True)

    depths = depths[times_no_clones]
    alleles = alleles[times_no_clones]

    freqs = estimate_frequencies(alleles, depths)

    for idx, mutation in enumerate(np.array(freqs, dtype=float)):
        masked_idx = np.nonzero(depths.iloc[idx, :].to_list())
        masked_times = times_no_clones[masked_idx]
        masked_mut = mutation[masked_idx]

        interpolating_function = create_interpolation_function(masked_times, masked_mut)
        fs = interpolating_function(times_no_clones)
        fs[fs < 0] = 0
        Ms += fs

    return (times_no_clones, Ms)


def average_trajectories(trajectories):
    avg_map = {}

    for ts, xs in trajectories:
        for t, x in zip(ts, xs):
            if t not in avg_map:
                avg_map[t] = {"x": 0.0, "n": 0.0}
            avg_map[t]["x"] += x
            avg_map[t]["n"] += 1

    for t in avg_map.keys():
        avg_map[t]["x"] /= avg_map[t]["n"]

    avg_ts = []
    avg_xs = []
    for t in sorted(avg_map.keys()):
        avg_ts.append(t)
        avg_xs.append(avg_map[t]["x"])

    return np.array(avg_ts), np.array(avg_xs)


def add_filter_to_timecourse(dfs: tuple, direction=None) -> tuple:
    codon_pos_dict = {}
    codon_direction_dict = {}

    for df in dfs:
        df.reset_index(drop=True, inplace=True)

    guide, depths, alleles = dfs

    guide["pairs"] = [f"{k}->{v}" for k, v in list(zip(guide["ref"], guide["alt"]))]
    guide["direction"] = [
        "AT->GC" if mut in at_to_gc else "GC->AT" for mut in guide["pairs"]
    ]

    df_codon_pos = guide.groupby("codon_pos")
    df_codon_and_direction = guide.groupby(["codon_pos", "direction"])

    if direction:
        # {0-AT->GC:(guide, depths, alleles), 1-AT->GC:(guide, depths, alleles), 2-AT->GC:(guide, depths, alleles)}
        for (codon_pos, direction), group in df_codon_and_direction:
            codon_direction_dict[f"{str(codon_pos)}-{direction}"] = (
                guide.loc[group.index],
                depths.loc[group.index],
                alleles.loc[group.index],
            )

        return codon_direction_dict

    else:  # only by codon_pos
        # {0:(guide, depths, alleles), 1:(guide, depths, alleles), 2:(guide, depths, alleles)}
        for codon_pos, group in df_codon_pos:
            codon_pos_dict[codon_pos] = (
                guide.loc[group.index],
                depths.loc[group.index],
                alleles.loc[group.index],
            )

        return codon_pos_dict


def estimate_frequencies(alts, depths):
    return alts * 1.0 / (depths + (depths == 0))


def create_interpolation_function(
    intpol_times, intpol_freqs, tmax=100000, kind="linear"
):
    padded_times = np.zeros(len(intpol_times) + 1)
    padded_freqs = np.zeros(len(intpol_times) + 1)
    padded_times[0 : len(intpol_times)] = intpol_times
    padded_freqs[0 : len(intpol_times)] = intpol_freqs
    padded_times[-1] = tmax
    padded_freqs[-1] = intpol_freqs[-1]

    interpolating_function = interp1d(
        padded_times, padded_freqs, kind=kind, bounds_error=True
    )

    return interpolating_function


def runner(files: list) -> dict:
    mutation_trajectories = {}
    avg_timecourse = {pos: [] for pos in POS_DIR_COLORS.keys()}

    for file in files:
        pop = str(file).split("/")[-1].split("_")[0]
        rs = format_timecourse(file)
        times, Ms = interpolate(rs)
        mutation_trajectories[pop] = (times, Ms)

        filtered_rs = add_filter_to_timecourse(rs, direction=True)
        if pop in MUTATORS:
            for pos, dfs in test.items():
                times, Ms = interpolate(dfs)
                avg_timecourse[pos].append((times, Ms))

    for pos, avg_dfs in avg_timecourse.items():
        avg_times, avg_Ms = average_trajectories(avg_dfs)
        # ax[1,1].plot(avg_times[:-7], avg_Ms[:-7], linestyle, color=pos_dir_colors[pos],markersize=1,markeredgewidth=0, linewidth=2.5)

    return [mutation_trajectories, avg_timecourse]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate mutations trajectories by codon"
    )
    parser.add_argument(
        "-p",
        "--path",
        help="Path to all timecourse files",
        type=str,
        metavar="PATH",
        required=True,
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
