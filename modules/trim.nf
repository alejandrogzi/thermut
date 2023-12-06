#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process TRIM {

    cpus 6

    input:
        tuple val(sample), path(fq)
        val metadata

    output:
        tuple val(sample), 
        path("*.paired.trim*.fastq.gz"), 
        emit: trim_paired

    script:
    """
    PREFIX=\$(awk -F'\t' '{print \$5}' <(grep ${sample} ${metadata}))
    TRIMDIR=\$(dirname \$(readlink -f \$(which trimmomatic)))/adapters

    if [[ \$PREFIX == CL* || \$PREFIX == batch* ]]; then
      ADAPTER="\$TRIMDIR/NexteraPE-PE.fa"
    elif [[ \$PREFIX == DL_LTE_rn1_* || \$PREFIX == ara* || \$PREFIX == SRR* || \$PREFIX == ERR* ]]; then
      ADAPTER="\$TRIMDIR/TruSeq3-PE-2.fa"
    else
      ADAPTER=''
    fi

    echo \$ADAPTER

    trimmomatic \\
    PE \\
    -threads ${task.cpus} \\
    ${fq[0]} ${fq[1]} \\
    ${sample}.paired.trim_1.fastq.gz ${sample}.unpaired.trim_1.fastq.gz \\
    ${sample}.paired.trim_2.fastq.gz ${sample}.unpaired.trim_2.fastq.gz \\
    -phred33 \\
    ILLUMINACLIP:\$ADAPTER:2:30:10:2:'false' LEADING:20 TRAILING:20

    rm ${sample}.unpaired.trim_1.fastq.gz ${sample}.unpaired.trim_2.fastq.gz
    """
}
