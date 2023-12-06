#!/usr/env/bin nextflow

nextflow.enable.dsl=2

process GET {
    input:
    val(run)

    output:
    tuple val(run), 
    path("*.fastq.gz"), 
    emit: fastq

    script:
    """
    ../../../bin/fastq-dl.py -a ${run} --cpus 1 --max-attempts 30 --sleep 120 
    """
}

workflow GET_FASTQ {
    take:
    runlist

    main:
    runs = Channel.fromPath(runlist).splitText().map{it.trim()}
    fastqs = GET(runs).fastq

    emit:
    fastqs
}
