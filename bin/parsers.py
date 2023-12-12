#!/usr/bin/env python3

# This is a script slightly modified from the original
# https://github.com/benjaminhgood/LTEE-metagenomic/blob/master/parse_file.py.
# It is used to parse the genome and annotation files.

# -----------------------------------------------------

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "alejandrogzi"
__version__ = "0.0.1"
__credits__ = ["Benjamin H Good"]

# -----------------------------------------------------

import numpy as np

# -----------------------------------------------------


GENOME = "../supp/REL606.gbk"
MASK = "../supp/REL606.L20.G15.P0.M35.mask.gd"


def parse_genome():
    with open(GENOME, "r") as file:
        for line in file:
            if line.startswith("ORIGIN"):
                break

        reference_sequences = (
            part.upper() for line in file for part in line.split() if part.isalpha()
        )
        reference_sequence = "".join(reference_sequences)

    return reference_sequence


def parse_mask_list():
    # Add masks calculated in Tenaillon et al (Nature, 2016)
    # Downloaded from barricklab/LTEE-Ecoli/reference/ github repository

    mask_start_positions = []
    mask_end_positions = []

    with open(MASK, "r") as file:
        file.readline()
        for line in file:
            items = line.split()
            start = int(items[4])
            length = int(items[5])
            mask_start_positions.append(start)
            mask_end_positions.append(start + length - 1)

    # Add masking of prophage elements (Methods section of Tenaillon et al (Nature, 2016))
    mask_start_positions.append(880528)
    mask_end_positions.append(904682)

    return np.array(mask_start_positions), np.array(mask_end_positions)


def parse_gene_list():
    features = set(["CDS", "gene", "tRNA", "rRNA", "repeat_region"])

    reference_sequence = parse_genome()
    observed_gene_names = set()

    gene_names = []
    start_positions = []
    end_positions = []
    promoter_start_positions = []
    promoter_end_positions = []
    gene_sequences = []
    strands = []

    with open(GENOME, "r") as file:
        line = file.readline()
        while line != "":
            items = line.split()
            feature = items[0]

            gene_name = ""
            feature_location = ""

            if feature == "CDS":
                feature_location = items[1]

                line = file.readline().strip()

                gene_name = ""
                locus_name = ""
                is_pseudo = False

                while line.split()[0] not in features:
                    if line.startswith("/gene"):
                        gene_name = line.split("=")[1].strip()[1:-1]
                    if line.startswith("/locus_tag"):
                        locus_name = line.split("=")[1].strip()[1:-1]
                    if line.startswith("/pseudo"):
                        is_pseudo = True

                    line = file.readline().strip()

                if gene_name == "":
                    gene_name = locus_name

                if is_pseudo:
                    gene_name = ""

                # done here

            elif feature == "tRNA" or feature == "rRNA":
                feature_location = items[1]
                # print feature_location

                while not line.strip().startswith("/gene"):
                    line = file.readline().strip()
                gene_name = line.split("=")[1].strip()[1:-1]
                gene_name = "%s:%s" % (feature, gene_name)

            else:
                # nothing to see here
                line = file.readline().strip()

            # If the element has a feature location string and a name
            # it should either be a gene, tRNA, or rRNA, so let's get details
            if feature_location != "" and gene_name != "":
                location_str = (
                    feature_location.lstrip("complement(").lstrip("join(").rstrip(")")
                )
                location_strs = location_str.split(",")

                for location_str in location_strs:
                    locations = [int(subitem) for subitem in location_str.split("..")]

                    gene_start = locations[0]
                    gene_end = locations[1]

                    if feature == "CDS":
                        gene_sequence = reference_sequence[gene_start - 1 : gene_end]
                    else:
                        gene_sequence = ""

                    strand = "forward"
                    promoter_start = (
                        gene_start - 100
                    )  # by arbitrary definition, we treat the 100bp upstream as promoters
                    promoter_end = gene_start - 1

                    if gene_sequence != "" and (not len(gene_sequence) % 3 == 0):
                        pass
                        print(gene_name, gene_start, "Not a multiple of 3")

                    if feature_location.startswith("complement"):
                        strand = "reverse"
                        promoter_start = gene_end + 1
                        promoter_end = gene_end + 100

                        if gene_sequence == "":
                            promoter_end = promoter_start

                    # record information

                    # first make sure gene name is unique
                    i = 1
                    old_gene_name = gene_name
                    while gene_name in observed_gene_names:
                        i += 1
                        gene_name = "%s_%d" % (old_gene_name, i)

                    start_positions.append(gene_start)
                    end_positions.append(gene_end)
                    promoter_start_positions.append(promoter_start)
                    promoter_end_positions.append(promoter_end)
                    gene_names.append(gene_name)
                    gene_sequences.append(gene_sequence)
                    strands.append(strand)
                    observed_gene_names.add(gene_name)

    # sort genes based on start position

    (
        gene_names,
        start_positions,
        end_positions,
        promoter_start_positions,
        promoter_end_positions,
        gene_sequences,
        strands,
    ) = (
        list(x)
        for x in zip(
            *sorted(
                zip(
                    gene_names,
                    start_positions,
                    end_positions,
                    promoter_start_positions,
                    promoter_end_positions,
                    gene_sequences,
                    strands,
                ),
                key=lambda pair: pair[1],
            )
        )
    )

    return (
        gene_names,
        np.array(start_positions),
        np.array(end_positions),
        np.array(promoter_start_positions),
        np.array(promoter_end_positions),
        gene_sequences,
        strands,
    )


def parse_repeat_list():
    repeat_names = []
    start_positions = []
    end_positions = []
    complements = []

    with open(GENOME, "r") as file:
        line = file.readline()
        while line != "":
            items = line.split()
            feature = items[0]

            if feature == "repeat_region":
                feature_location = items[1]

                # Get name of mobile element
                repeat_name = "unknown"

                line = file.readline()
                while line.strip()[0] == "/":
                    if line.strip().startswith("/mobile_element"):
                        repeat_name = line.split("=")[1].strip()[1:-1]
                    line = file.readline()

                # Finished at next non '/' entry, presumably next feature

                if feature_location.startswith("complement"):
                    complement = True
                else:
                    complement = False

                location_str = (
                    feature_location.lstrip("complement(").lstrip("join(").rstrip(")")
                )
                location_strs = location_str.split(",")
                for location_str in location_strs:
                    locations = [int(subitem) for subitem in location_str.split("..")]
                    start_positions.append(locations[0])
                    end_positions.append(locations[1])
                    repeat_names.append(repeat_name)
                    complements.append(complement)

            else:
                line = file.readline()

    return repeat_names, np.array(start_positions), np.array(end_positions), complements
