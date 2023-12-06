#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process GROUP {
    input:
    val sample
    val dir

    output:
    tuple val(sample), val("${dir}/${sample}_1.fastq.gz"), val("${dir}/${sample}_2.fastq.gz"), emit: fastq

    script:
    """
    """
}


workflow {
    group = params.group // file containing list of runs
    dir = params.dir // directory containing fastq files
    runs = Channel.fromPath(group).splitText().map{it.trim()}
    fastqs = GROUP(runs, dir).fastq

    fastqs.view()
}
