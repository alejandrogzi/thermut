#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process BRESEQ {

  cpus 4

  input:
  tuple val(sample), path(fqs)
  path(genome)

  output:
  path('data/*.bam'), emit: bam

  script:
  """
  breseq -p -j ${task.cpus} --brief-html-output \\
  --polymorphism-reject-indel-homopolymer-length 0 \\
  --polymorphism-reject-surrounding-homopolymer-length 0 \\
  --polymorphism-score-cutoff 2 \\
  -r ${genome} ${fqs}

  mv data/reference.bam data/${sample}.bam
  """
}
