#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process ANNOTATE {

  input:
  tuple val(group), path(likelihood)

  output:
  tuple val(group), path('*.txt'), emit annotate_ch

  script:
  """
  ../../../bin/annotate.py --likelihood ${likelihood} --population ${group} --out ${group}.annotated.timecourse.txt
  """
}
