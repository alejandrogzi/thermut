#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// need to check if snp.nf outputs files in order.
// the channel should be ordered like:
// depth.txt -> snp.txt 
// if this is not true, the C++ binary will make
// things wrong.
process PVAL {

  input:
  tuple val(group), path(mutations)

  output:
  tuple val(group), path('*.txt'), emit: sign_mutations

  script:
  """
  g++ --std=c++11 -O3 ../../../pv/main.cpp -o ../../../pv/calculate_pval
  wait

  cat ${mutations} | ../../../pv/calculate_pval >  ${group}_likelihood_timecourse.txt
  """
}
