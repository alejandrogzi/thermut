#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process GROUP {
    input:
    val sample
    val dir

    output:
    tuple val(sample), 
    val("${dir}/${sample}_1.fastq.gz"), 
    val("${dir}/${sample}_2.fastq.gz"), 
    emit: fastq

    script:
    """
    """
}

workflow MAKE_GROUP {
    take:
    group // file with list of samples
    dir // directory with fastq files

    main:
    runs = Channel.fromPath(group).splitText().map{it.trim()}
    fastqs = GROUP(runs, dir).fastq.map{tuple(it[0], it[1..-1])}

    emit:
    fastqs
}
