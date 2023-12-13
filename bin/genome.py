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

import parsers as gp

# -----------------------------------------------------

GENOME_LENGTH = 4629812
NUCLEOTIDE = {"A": "T", "T": "A", "G": "C", "C": "G"}
CODON = {
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "AAT": "N",
    "AAC": "N",
    "GAT": "D",
    "GAC": "D",
    "TGT": "C",
    "TGC": "D",
    "CAA": "Q",
    "CAG": "Q",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "CAT": "H",
    "CAC": "H",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "AAA": "K",
    "AAG": "K",
    "ATG": "M",
    "TTT": "F",
    "TTC": "F",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "TGG": "W",
    "TAT": "Y",
    "TAC": "Y",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TAA": "!",
    "TGA": "!",
    "TAG": "!",
}


def calculate_substitutions() -> (dict, dict, dict, list):
    # calculate number of synonymous opportunities for each codon
    codon_synonymous_opportunity_table = {}
    for codon in CODON.keys():
        codon_synonymous_opportunity_table[codon] = {}
        for i in range(0, 3):
            codon_synonymous_opportunity_table[codon][
                i
            ] = (
                -1
            )  # since G->G is by definition synonymous, but we don't want to count it
            codon_list = list(codon)
            for base in ["A", "C", "T", "G"]:
                codon_list[i] = base
                new_codon = "".join(codon_list)
                if CODON[codon] == CODON[new_codon]:
                    # synonymous!
                    codon_synonymous_opportunity_table[codon][i] += 1

    bases = set(["A", "C", "T", "G"])
    substitutions = []
    for b1 in bases:
        for b2 in bases:
            if b2 == b1:
                continue

            substitutions.append("%s->%s" % (b1, b2))

    codon_synonymous_substitution_table = {}
    codon_nonsynonymous_substitution_table = {}
    for codon in CODON.keys():
        codon_synonymous_substitution_table[codon] = [[], [], []]
        codon_nonsynonymous_substitution_table[codon] = [[], [], []]

        for i in range(0, 3):
            reference_base = codon[i]

            codon_list = list(codon)
            for derived_base in ["A", "C", "T", "G"]:
                if derived_base == reference_base:
                    continue
                substitution = "%s->%s" % (reference_base, derived_base)
                codon_list[i] = derived_base
                new_codon = "".join(codon_list)
                if CODON[codon] == CODON[new_codon]:
                    # synonymous!
                    codon_synonymous_substitution_table[codon][i].append(substitution)
                else:
                    codon_nonsynonymous_substitution_table[codon][i].append(
                        substitution
                    )

    return (
        codon_synonymous_opportunity_table,
        codon_synonymous_substitution_table,
        codon_nonsynonymous_substitution_table,
        substitutions,
    )


def calculate_reverse_complement_sequence(seq: str) -> str:
    return "".join(NUCLEOTIDE[base] for base in seq[::-1])


def annotation_map() -> (dict, dict, dict):
    (
        codon_synonymous_opportunity_table,
        codon_synonymous_substitution_table,
        codon_nonsynonymous_substitution_table,
        substitutions,
    ) = calculate_substitutions()

    (
        gene_names,
        gene_start_positions,
        gene_end_positions,
        promoter_start_positions,
        promoter_end_positions,
        gene_sequences,
        strands,
    ) = gp.parse_gene_list()
    (
        repeat_names,
        repeat_start_positions,
        repeat_end_positions,
        repeat_complements,
    ) = gp.parse_repeat_list()
    mask_start_positions, mask_end_positions = gp.parse_mask_list()

    position_gene_map = {}
    gene_position_map = {}

    num_masked_sites = 0

    # first mark things that are repeats
    # this takes precedence over all other annotations
    for start, end in zip(repeat_start_positions, repeat_end_positions):
        for position in range(start, end + 1):
            if position not in position_gene_map:
                position_gene_map[position] = "repeat"
                num_masked_sites += 1

    # then mark masked things
    for start, end in zip(mask_start_positions, mask_end_positions):
        for position in range(start, end + 1):
            if position not in position_gene_map:
                position_gene_map[position] = "repeat"
                num_masked_sites += 1

    # then greedily annotate genes at remaining sites
    for gene_name, start, end in zip(
        gene_names, gene_start_positions, gene_end_positions
    ):
        for position in range(start, end + 1):
            if position not in position_gene_map:
                position_gene_map[position] = gene_name
                if gene_name not in gene_position_map:
                    gene_position_map[gene_name] = []
                gene_position_map[gene_name].append(position)

    # remove 'partial' genes that have < 10bp unmasked sites
    for gene_name in list(sorted(gene_position_map.keys())):
        if len(gene_position_map[gene_name]) < 10:
            for position in gene_position_map[gene_name]:
                position_gene_map[position] = "repeat"
            del gene_position_map[gene_name]

    # count up number of synonymous opportunities
    effective_gene_synonymous_sites = {}
    effective_gene_nonsynonymous_sites = {}

    substitution_specific_synonymous_sites = {
        substitution: 0 for substitution in substitutions
    }
    substitution_specific_nonsynonymous_sites = {
        substitution: 0 for substitution in substitutions
    }

    for gene_name, start, end, gene_sequence, strand in zip(
        gene_names, gene_start_positions, gene_end_positions, gene_sequences, strands
    ):
        if gene_name not in gene_position_map:
            continue

        if strand == "forward":
            oriented_gene_sequence = gene_sequence
        else:
            oriented_gene_sequence = calculate_reverse_complement_sequence(
                gene_sequence
            )

        for position in gene_position_map[gene_name]:
            if gene_name not in effective_gene_synonymous_sites:
                effective_gene_synonymous_sites[gene_name] = 0
                effective_gene_nonsynonymous_sites[gene_name] = 0

            if gene_name.startswith("tRNA") or gene_name.startswith("rRNA"):
                pass

            else:
                # calculate position in gene
                if strand == "forward":
                    position_in_gene = position - start
                else:
                    position_in_gene = end - position

                # calculate codon start
                codon_start = int(position_in_gene / 3) * 3
                codon = gene_sequence[codon_start : codon_start + 3]
                position_in_codon = position_in_gene % 3

                # print gene_name, start, end, position, codon,position_in_codon

                effective_gene_synonymous_sites[gene_name] += (
                    codon_synonymous_opportunity_table[codon][position_in_codon] / 3.0
                )
                effective_gene_nonsynonymous_sites[gene_name] += (
                    1
                    - codon_synonymous_opportunity_table[codon][position_in_codon] / 3.0
                )

                for substitution in codon_synonymous_substitution_table[codon][
                    position_in_codon
                ]:
                    substitution_specific_synonymous_sites[substitution] += 1

                for substitution in codon_nonsynonymous_substitution_table[codon][
                    position_in_codon
                ]:
                    substitution_specific_nonsynonymous_sites[substitution] += 1

    substitution_specific_synonymous_fraction = {
        substitution: substitution_specific_synonymous_sites[substitution]
        * 1.0
        / (
            substitution_specific_synonymous_sites[substitution]
            + substitution_specific_nonsynonymous_sites[substitution]
        )
        for substitution in substitution_specific_synonymous_sites.keys()
    }

    # then annotate promoter regions at remaining sites
    for gene_name, start, end in zip(
        gene_names, promoter_start_positions, promoter_end_positions
    ):
        for position in range(start, end + 1):
            if position not in position_gene_map:
                # position hasn't been annotated yet

                if gene_name not in gene_position_map:
                    # the gene itself has not been annotated
                    # so don't annotate the promoter
                    pass
                else:
                    position_gene_map[position] = gene_name
                    gene_position_map[gene_name].append(position)

    # calculate effective gene lengths
    effective_gene_lengths = {
        gene_name: len(gene_position_map[gene_name])
        - effective_gene_synonymous_sites[gene_name]
        for gene_name in gene_position_map.keys()
    }
    effective_gene_lengths["synonymous"] = sum(
        [
            effective_gene_synonymous_sites[gene_name]
            for gene_name in gene_position_map.keys()
        ]
    )
    effective_gene_lengths["nonsynonymous"] = sum(
        [
            effective_gene_nonsynonymous_sites[gene_name]
            for gene_name in gene_position_map.keys()
        ]
    )
    effective_gene_lengths["masked"] = num_masked_sites
    effective_gene_lengths["noncoding"] = (
        GENOME_LENGTH
        - effective_gene_lengths["synonymous"]
        - effective_gene_lengths["nonsynonymous"]
        - effective_gene_lengths["masked"]
    )

    return (
        position_gene_map,
        effective_gene_lengths,
        substitution_specific_synonymous_fraction,
    )


def is_repeat(pos, pos_map):
    if pos in pos_map:
        if pos_map[pos] == "repeat":
            return True
    return False


def annotate_gene(position, position_gene_map):
    if position in position_gene_map:
        gene_name = position_gene_map[position]
    else:
        gene_name = "intergenic"

    return gene_name


def annotate_variant(position, allele, gene_data, position_gene_map):
    (
        gene_names,
        gene_start_positions,
        gene_end_positions,
        promoter_start_positions,
        promoter_end_positions,
        gene_sequences,
        strands,
    ) = gene_data

    # get gene
    gene_name = annotate_gene(position, position_gene_map)

    if allele.startswith("Depth"):
        var_type = "unknown"
    elif allele[1:3] == "->":
        # a SNP, so annotate it
        if gene_name == "intergenic":
            var_type = "noncoding"
        elif gene_name == "repeat":
            var_type = "repeat"
        else:
            # must be in a real gene
            # so get it
            i = gene_names.index(gene_name)

            gene_start_position = gene_start_positions[i]
            gene_end_position = gene_end_positions[i]
            promoter_start_position = promoter_start_positions[i]
            promoter_end_position = promoter_end_positions[i]
            gene_sequence = gene_sequences[i]
            strand = strands[i]

            if position < gene_start_position or position > gene_end_position:
                # var_type='promoter'
                var_type = "noncoding"  # (promoter)
            else:
                if gene_name.startswith("tRNA") or gene_name.startswith("rRNA"):
                    var_type = "noncoding"
                else:
                    # calculate position in gene
                    if strand == "forward":
                        position_in_gene = position - gene_start_position
                        oriented_gene_sequence = gene_sequence
                        new_base = allele[3]
                    else:
                        position_in_gene = gene_end_position - position
                        oriented_gene_sequence = calculate_reverse_complement_sequence(
                            gene_sequence
                        )
                        new_base = NUCLEOTIDE[allele[3]]

                    # calculate codon start
                    codon_start = int(position_in_gene / 3) * 3
                    codon = oriented_gene_sequence[codon_start : codon_start + 3]
                    codon_list = list(codon)
                    position_in_codon = position_in_gene % 3
                    codon_list[position_in_codon] = new_base
                    new_codon = "".join(codon_list)
                    if CODON[codon] == CODON[new_codon]:
                        var_type = "synonymous"
                    else:
                        if CODON[new_codon] == "!":
                            var_type = "nonsense"
                        else:
                            var_type = "missense"
    else:
        var_type = "unknown"

    return gene_name, var_type
