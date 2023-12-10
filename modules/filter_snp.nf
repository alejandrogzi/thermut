#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process FILTER_SNP {

  input:
  tuple val(group), path(timecourse)

  output:
  tuple val(group), path('*.txt'), emit: snp_ch

  script:
  """
  ../../../bin/filter_snp.py -i ${timecourse} -p ${group}
  """
}
