#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TIMECOURSE {

  cpus 10

  input:
  path(bams)
  path(fasta)
  val group

  output:
  tuple val(group), path("*.timecourse"), emit: timecourse_ch

  script:
  def name = group.split("/")[-1][0..1]
  """
  for bam in ${bams} ; do
    samtools index \${bam}
  done

  parallel \\
  --keep-order \\
  --jobs ${task.cpus} \\
  --colsep '\t' \\
  'bcftools mpileup -Ou -d 1000 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -q10 -f ${fasta} -r {1}:{2}-{3} ${bams} | bcftools call --ploidy 1 -mv -Oz -o output.{1}.{2}.{3}.vcf.gz' :::: ../../../supp/intervals.bed

  bcftools concat -Oz -o output.vcf.gz output*.vcf.gz
  rm output.REL*.vcf.gz

  ../../../bin/timecourse.py --vcf output.vcf.gz --output ${name}.timecourse
  """
}

