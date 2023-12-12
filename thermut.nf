#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Usage:
// nextflow run main.nf --dir /path/to/dir --meta /path/to/meta.csv --group /path/to/group.csv --gbk /path/to/gbk --fasta /path/to/fasta
//
// Where:
// dir = directory containing fastq files
// meta = metadata file (./meta/metadata.txt)
// group = group file (a list of samples per population produced by ./grouper.sh)
// gbk = reference genbank file (./supp/REL606.6.gbk)
// fasta = reference fasta file (./supp/REL606.6.fasta)


// modules
include { TRIM } from './modules/trim'
include { BRESEQ } from './modules/breseq'
include { TIMECOURSE } from './modules/timecourse'
include { FILTER_SNP } from './modules/filter_snp'
include { PVAL } from './modules/pval'
include { ANNOTATE } from './modules/annotate'

// subworkflows
include { MAKE_GROUP } from './subworkflows/group'
include { GET } from './subworkflows/getter'


workflow {
  dir = params.dir
  meta = params.meta
  group = params.group
  gbk = params.gbk
  fasta = params.fasta

  fqs = MAKE_GROUP(group, dir).fastqs
  trim_fqs = TRIM(fqs, meta).trim_paired
  bams = BRESEQ(trim_fqs, gbk).bam.collect()
  timecourse = TIMECOURSE(bams, fasta, group)
  snp = FILTER_SNP(timecourse)
  pval = PVAL(snp)
  annotate = ANNOTATE(pval)
}
