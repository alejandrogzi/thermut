#!/usr/bin/env python3

# This script takes a timecourse file and
# filters out SNPs that are in repeat regions
#
# Distinct SNPs at the same site are treated as two separate mutations
# and INDELS are explicitly ignored.
#
# Usage:
# python3 annotate.py --likelihood <file.likelihood> --population <population>

# -----------------------------------------------------

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "alejandrogzi"
__version__ = "0.0.1"
__credits__ = ["Benjamin H Good"]

# -----------------------------------------------------

import numpy as np
import argparse
import pandas as pd
import genome
import parsers
import utils

# -----------------------------------------------------

araA = 70867  # , araA, T->C @ 70867
recD = 2847052  # , recD, A->G @ 2847052
FDR = 0.05
MUTATORS = ["m1", "m2", "m3", "m4", "p3", "p6"]
NONMUTATORS = ["m5", "m6", "p1", "p2", "p4", "p5"]
POS_MAP, GENE_LENGTHS, SUBS_FRACTION = genome.annotation_map()
GENES = parsers.parse_gene_list()
REPEATS = parsers.parse_repeat_list()
MASKS = parsers.parse_mask_list()


def main():
    args = parse_args()
    annotate(args.likelihood, args.population)


def set_pval(likelihood: str, population: str) -> float:
    """Set p-value for a given group of populations [mutators, nonmutators]

    @param: likelihood
    @param: population
    @return: p-value
    """

    pvas = np.array([])
    with open(likelihood, "r") as file:
        for line in file:
            fields = line.strip().split(",")
            pos = int(fields[1])
            pva = float(fields[6].split()[1])  # combined p-value
            corr = float(fields[6].split()[0])

            if pos == 0:
                continue

            if genome.annotate_gene(pos, POS_MAP) == "repeat":
                continue

            np.append(pvas, pva)

    for v in sorted(pvas, reverse=True):
        if len(pvas) * v / (pvas <= v).sum() <= FDR:
            return v


def annotate(likelihood: str, population: str, out: str) -> None:
    """Annotate SNPs in a given group of populations [mutators, nonmutators]

    @param: likelihood
    @param: population
    """

    THRESHOLD_PVAL = set_pval(likelihood, population)
    num_total = 0  # trajectories
    num_passed = 0  # trajectories passed

    annotated_mutations = []
    header = None
    printed_header = False
    loaded_avg_depths = False

    with open(likelihood, "r") as file:
        for line in file:
            fields = line.strip().split(",")
            chr = fields[0].strip()
            pos = int(fields[1])
            allele = fields[2].strip()
            times = np.array(fields[3].strip().split(), dtype=float) * 1000
            alts = np.array(fields[4].strip().split(), dtype=float)
            depths = np.array(fields[5].strip().split(), dtype=float)
            stat = float(fields[6][0])
            pval = float(fields[6][1])

            deletion_idx, fold_reduction, deletion_pvaue = tuple(
                [float(subitem) for subitem in fields[7].split()]
            )
            deletion_idx = int(deletion_idx)

            duplication_idx, fold_increase, duplication_pvaue = tuple(
                [float(subitem) for subitem in fields[8].split()]
            )
            duplication_idx = int(duplication_idx)

            of_times = times[times < 1000000]
            of_alts = alts[times < 1000000]
            of_depths = depths[times < 1000000]

            clone_times = times[times > 1000000] - 1000000
            clone_alts = alts[times > 1000000]
            clone_depths = depths[times > 1000000]

            # load average depths if not leaded yet
            if not loaded_avg_depths:
                pop_avg_depths = depths
                clone_avg_depths = clone_depths
                loaded_avg_depths = True

            # create header if not created yet
            if header == None:
                print_strings = [
                    "Position",
                    "Gene",
                    "Allele",
                    "Annotation",
                    "Test statistic",
                    "P-value",
                    "Deletion index",
                    "Fold reduction",
                    "Deletion P-value",
                    "Duplication index",
                    "Fold increase",
                    "Duplication pvaue",
                    "Passed?",
                ]
                for t in zip(of_times):
                    print_strings.append("AC:%d" % t)
                    print_strings.append("DP:%d" % t)

                header = ",".join(print_strings)

            # annotate mutation
            gene_name, var_type = genome.annotate_variant(pos, allele, GENES, POS_MAP)

            # if pvaue is lower than threshold and not a weird insertion gene
            passed_str = "FAIL"

            if (
                (pval <= THRESHOLD_PVAL)
                and gene_name != "repeat"
                and (pos != araA)
                and (pos != recD)
                and (alts[0] * 1.0 / (depths[0] + (depths[0] == 0)) < 0.2)
            ):
                # would otherwise pass
                # check if clone fail

                # determine whether clone data suggests that the mutation
                # is a duplication. do not pass these

                # first estimate frequencies at good timepoints
                good_idxs, filtered_alts, filtered_depths = utils.mask_timepoints(
                    times,
                    alts,
                    depths,
                    var_type,
                    deletion_idx,
                    fold_reduction,
                    deletion_pvaue,
                )
                freqs = utils.estimate_frequencies(filtered_alts, filtered_depths)
                masked_times = times[good_idxs]
                masked_freqs = freqs[good_idxs]

                masked_depth_ratios = depths[good_idxs] / pop_avg_depths[good_idxs]

                interpolation_function = utils.create_interpolation_function(
                    masked_times, masked_freqs, tmax=100000
                )

                (
                    masked_clone_times,
                    masked_clone_freqs,
                ) = utils.estimate_clone_frequencies(
                    clone_times, clone_alts, clone_depths
                )

                if len(masked_clone_times) == 0:
                    # clone fail.. there are no points with actual reads.. this is an issue
                    clone_fail = True

                else:
                    # we have some clone timepoints to work with
                    clone_fail = False

                    # try to detect duplication without clone information
                    if (masked_freqs > 0.25).sum() > 10:  # sticks around for a while
                        if (
                            masked_depth_ratios[masked_freqs > 0.25].mean()
                            / (masked_depth_ratios[:10].mean())
                            > 1.5
                        ):  # evidence of duplication
                            if (
                                (masked_freqs > 0.25) * (masked_freqs < 0.75)
                            ).sum() * 1.0 / (
                                masked_freqs > 0.25
                            ).sum() > 0.9:  # never really fixes
                                if (
                                    masked_clone_freqs < 0.9
                                ).all():  # just to make sure we aren't picking a trajectory that has fixed in a clone
                                    clone_fail = True

                    # now try to detect deletion with clone information

                    # frequency in population at time the clone is called
                    pop_freqs = interpolation_function(masked_clone_times)

                    # number of timepoints where mutation is at intermediate frequency in a clone
                    intermediate_clone_idxs = (masked_clone_freqs > 0.10) * (
                        masked_clone_freqs < 0.7
                    )
                    num_intermediate_clones = intermediate_clone_idxs.sum()
                    if intermediate_clone_idxs.sum() >= 4:
                        if (
                            (pop_freqs > 0.2)
                            * (pop_freqs < 0.7)
                            * intermediate_clone_idxs
                        ).sum() * 1.0 / intermediate_clone_idxs.sum() >= 0.5:
                            clone_fail = True

                    if ((masked_clone_freqs < 0.6) * (pop_freqs < 0.6)).all() and (
                        (masked_clone_freqs > 0.1) * (pop_freqs > 0.1)
                    ).any():
                        clone_fail = True

                    # see if there is evidence for a duplication
                    # calculate clone depth changes specifically where coverage is ok
                    clone_freqs = (
                        clone_alts * 1.0 / (clone_depths + (clone_depths == 0))
                    )
                    clone_depth_fold_changes = utils.estimate_depth_fold_changes(
                        clone_avg_depths, clone_depths
                    )
                    clone_idxs = np.array(
                        [t in masked_clone_times for t in clone_times]
                    )
                    clone_depth_fold_changes = clone_depth_fold_changes[clone_idxs]
                    num_duplicated_intermediates = (
                        clone_depth_fold_changes[
                            (masked_clone_freqs > 0.4) * (masked_clone_freqs < 0.8)
                        ]
                        > 0.4
                    ).sum() * 1.0

                    if (
                        (num_duplicated_intermediates > 0)
                        and (masked_clone_freqs < 0.9).all()
                        and (masked_freqs < 0.9).all()
                    ):
                        clone_fail = True

                    if masked_freqs[0] > 0.1:
                        clone_fail = True

                if not clone_fail:
                    passed_str = "PASS"
                    num_passed += 1

            # print to CSV file
            print_strings = [
                str(pos),
                gene_name,
                allele,
                var_type,
                str(stat),
                str(pval),
                str(deletion_idx),
                str(fold_reduction),
                str(deletion_pvaue),
                str(duplication_idx),
                str(fold_increase),
                str(duplication_pvaue),
                passed_str,
            ]
            for alt, depth in zip(alts, depths):
                print_strings.append(str(alt))
                print_strings.append(str(depth))

            annotated_mutations.append((pos, ",".join(print_strings)))

    # sort mutations by Position
    annotated_mutations.sort(key=lambda x: x[0])

    # write to file
    with open(out, "w") as file:
        file.write(header + "\n")
        for pos, line in annotated_mutations:
            file.write(line + "\n")

    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate SNPs in a given group of populations [mutators, nonmutators]"
    )
    parser.add_argument(
        "-l", "--likelihood", help="Likelihood file", required=True, type=str
    )
    parser.add_argument(
        "-p", "--population", help="Population", required=True, type=str
    )
    parser.add_argument("-o", "--output", help="Output file", required=True, type=str)
    return parser.parse_args()


if __name__ == "__main__":
    main()
