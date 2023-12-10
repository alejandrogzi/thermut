#!/usr/bin/env python3

# This script takes a timecourse file and
# filters out SNPs that are in repeat regions
#
# Distinct SNPs at the same site are treated as two separate mutations
# and INDELS are explicitly ignored.
#
# Usage:
# python3 filter_snp.py --timecourse <file.timecourse> --population <population>

# -----------------------------------------------------

__author__ = 'Alejandro Gonzales-Irribarren'
__email__ = 'jose.gonzalesdezavala1@unmsm.edu.pe'
__github__ = 'alejandrogzi'
__version__ = '0.0.1'
__credits__ = ['Benjamin H Good']

# -----------------------------------------------------


import lensky.pipeline.bin.genome as genome
import numpy as np
import argparse

# -----------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--timecourse", required=True ,help="Timecourse file")
    parser.add_argument("-p", "--population", required=True, help="Name of population")
    args = parser.parse_args()
    return args

def filter_snp(timecourse: str, pop: str):
    f_snp = open(f"./snp_{pop}.txt", "w")
    f_depth = open(f"./depth_{pop}.txt", "w")
    
    re_depth = []
    re_time = []

    pos_map, gene_lengths, substitution_fraction = genome.annotation_map()

    with open(timecourse, "r") as file:
        for line in file:
            fields = line.split(",")
            pos = int(fields[1])
            allele = fields[2].strip()
            times = fields[3]

            if genome.is_repeat(pos, pos_map):
                continue

            f_snp.write(line)

            times = np.fromstring(fields[3], sep=" ")
            depths = np.fromstring(fields[5], sep=" ")

            re_depth.append(depths)
            re_time.append(times)

        re_depth = np.vstack(re_depth)
        median_depths = np.median(re_depth, axis=0)
        alts = [0]*len(median_depths)

        depth_line = f'REL606,0,Depth,{" ".join(map(str, re_time[0]))},{" ".join(map(str, alts))},{" ".join(map(str, median_depths))}\n'
        f_depth.write(depth_line)


def main():
    args = parse_args()
    filter_snp(args.timecourse, args.population)


if __name__ == "__main__":
    main()