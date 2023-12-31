#!/usr/bin/env python3

# This script takes an vcf.gz and calculates a timecourse for
# mutations called. This is a re-implementation of the original:
# https://github.com/benjaminhgood/LTEE-metagenomic/blob/master/cluster_scripts/create_timecourse.py
#
# Distinct SNPs at the same site are treated as two separate mutations
# and INDELS are explicitly ignored.
#
# Usage:
# python3 timecourse.py --vcf <vcf.gz> --out <output>

# -----------------------------------------------------

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "alejandrogzi"
__version__ = "0.0.1"
__credits__ = ["Benjamin H Good"]

# -----------------------------------------------------

import numpy as np
import pandas as pd
import argparse
import gzip

# -----------------------------------------------------

FIRST_POSITION = 0
LAST_POSITION = 4629812
METADATA = "../../../meta/metadata.txt"
MIN_ALLELE_FREQ = 2
MIN_DEPTH = 10
MIN_EMP_FREQ = 0.05


def main():
    args = parse_args()
    vcf = args.vcf
    output = args.output
    timecourse(vcf, output)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Calculate timecourse for mutations at each site in the genome"
    )
    parser.add_argument(
        "--vcf", required=True, type=str, metavar="<vcf.gz>", help="vcf.gz file"
    )
    parser.add_argument(
        "--output", required=True, type=str, metavar="<output>", help="output file"
    )
    return parser.parse_args()


def timecourse(vcf: str, output: str) -> None:
    vcf = gzip.open(vcf, "r")
    meta = pd.read_csv(METADATA, sep="\t").set_index("Run")
    out = open(output, "w")

    for line in vcf:
        if line.startswith(b"##"):
            continue
        elif line.startswith(b"#CHROM"):
            headers = line.split(b"\t")
            continue

        fields = line.split(b"\t")

        if fields[7].startswith(b"INDEL"):
            continue

        chr = fields[0].decode("utf-8")  # chromosome
        pos = int(fields[1])  # genome position
        ref = fields[3]  # reference allele
        alt = fields[4]  # alternate allele
        sps = len(fields[9:])  # number of samples

        alleles = {ref: [0] * sps, alt: [0] * sps}
        depth = []
        times = []

        if pos < FIRST_POSITION:
            continue
        elif pos > LAST_POSITION:
            break

        for idx in range(9, len(fields)):
            id = headers[idx].split(b".")[0].decode("utf-8")
            times.append(meta.loc[id, "time"])
            depth.append(
                int(fields[idx].split(b":")[2].decode("utf-8"))
            )  # registers read depths for each sample
            gt = fields[idx].split(b":")[0]  # registers a sample genotype

            support = fields[idx].split(b":")[-1].split(b",")
            ref_support = support[0]
            alt_support = support[1]

            alleles[ref][idx - 9] += int(ref_support.decode("utf-8"))
            alleles[alt][idx - 9] += int(alt_support.decode("utf-8"))

        depth_map = {}

        for idx in range(0, len(depth)):
            if times[idx] not in depth_map:
                depth_map[times[idx]] = 0
            depth_map[times[idx]] += depth[idx]

        merged_times = np.array([t for t in sorted(depth_map.keys())])
        merged_depths = np.array([depth_map[t] for t in merged_times])
        alt_map = {}

        for allele in alleles.keys():
            if allele != ref:
                if allele not in alt_map:
                    alt_map[allele] = {t: 0 for t in depth_map.keys()}

                for i in range(0, len(depth)):
                    alt_map[allele][times[i]] += alleles[allele][i]

        merged_alts = {}
        for allele_key in alt_map.keys():
            merged_alts[allele_key] = np.array(
                [alt_map[allele_key][t] for t in merged_times]
            )

        for allele_key in merged_alts.keys():
            mutation = "%s->%s" % (ref.decode("utf-8"), allele_key.decode("utf-8"))
            allele_depths = merged_depths
            allele_alts = merged_alts[allele_key]

        allele_freq = (allele_alts >= 2).sum() >= MIN_ALLELE_FREQ
        sample_threshold = (
            (allele_alts >= MIN_ALLELE_FREQ)
            * (allele_depths >= MIN_DEPTH)
            * (allele_alts >= MIN_EMP_FREQ * allele_depths)
        ).sum() > 0

        if allele_freq and sample_threshold:
            outline = ",".join(
                [
                    chr,
                    str(pos),
                    mutation,
                    " ".join(str(t) for t in merged_times),
                    " ".join(str(a) for a in allele_alts),
                    " ".join(str(d) for d in allele_depths),
                ]
            )
            out.write(outline + "\n")


if __name__ == "__main__":
    main()
