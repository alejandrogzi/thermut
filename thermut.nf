#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// modules
include { TRIM } from './modules/trim'
include { BRESEQ } from './modules/breseq'
include { TIMECOURSE } from './modules/timecourse'
include { FILTER_SNP } from './modules/snp'

// subworkflows
include { MAKE_GROUP } from './subworkflows/group'
include { GET } from './subworkflows/getter'


workflow {
  dir = params.dir
  meta = params.meta
  group = params.group
  gbk = params.gbk
  //fasta = params.fasta

  fqs = MAKE_GROUP(group, dir).fastqs
  trim_fqs = TRIM(fqs, meta).trim_paired
  bams = BRESEQ(trim_fqs, gbk).bam.collect()
 // timecourse = TIMECOURSE(bams, fasta, group)
 // snp = FILTER_SNP(timecourse)
}
