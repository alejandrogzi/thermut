#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TIMECOURSE {

  cpus 10

  input:
  path(bams)
  path(fasta)

  output:

  script:
  """
  parallel \\
  --jobs ${task.cpus} \\
  --colsep '\t' \\
  'bcftools mpileup -Ou -d 1000 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -q10 -f ${fasta} -r {1}:{2}-{3} ${bams} | bcftools call --ploidy 1 -mv -Oz -o output.{1}.{2}.{3}.vcf.gz' :::: intervals.bed

  bcftools concat -Oz -o output.vcf.gz output*.vcf.gz
  rm output.REL*.vcf.gz

  ./timecourse.py --vcf ${sample}.mpileup --out ${sample}.timecourse 
  """
}




