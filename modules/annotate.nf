#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process ANNOTATE {

  input:
  tuple val(group), path(likelihood)

  output:
  tuple val(group), path('*.txt'), emit: annotate_ch

  script:
  def name = group.split('/')[-1][0..1]
  """
  ../../../bin/annotate.py --likelihood ${likelihood} --population ${name} --out ${name}.annotated.timecourse.txt
  """
}
