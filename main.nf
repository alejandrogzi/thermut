#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process INDEX {
    input:
    file(fasta)

    output:
    path(bwa), emit: idx

    script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir bwa
    bwa \\
        index \\
        -p bwa/${prefix} \\
        ${fasta}
    """
}


process TRIMMOMATIC {

    cpus 6

    input:
        tuple val(sample), path(fastq)
        val dir

    output:
        tuple val(sample), path("*.paired.trim*.fastq.gz"), emit: trim_paired

    script:
    """
    trimmomatic \\
    PE \\
    -threads ${task.cpus} \\
    ${fastq[0]} ${fastq[1]} \\
    ${sample}.paired.trim_1.fastq.gz ${sample}.unpaired.trim_1.fastq.gz \\
    ${sample}.paired.trim_2.fastq.gz ${sample}.unpaired.trim_2.fastq.gz \\
    -phred33 \\
    SLIDINGWINDOW:4:20 MINLEN:30

    rm ${sample}.unpaired.trim_1.fastq.gz ${sample}.unpaired.trim_2.fastq.gz
    """
}

process DEDUP {

    cpus 2

    input:
    path(idx)
    tuple val(sample), file(reads)

    output:
    tuple val(sample), path("${sample}.dedup.bam"), emit: dedup_bam

    script:
    """
    INDEX=\$(find -L ./ -name "*.amb" | sed 's/\\.amb\$//')
    
    bwa mem \$INDEX ${reads} -t ${task.cpus} | \\
    samtools view -b | samtools sort --threads 5 -o ${sample}.bam \\
    && picard MarkDuplicates -I ${sample}.bam -O ${sample}.dedup.bam \\
    -M metrics.txt --REMOVE_DUPLICATES true && \\
    samtools index ${sample}.dedup.bam
    
    rm ${sample}.bam
    """
}

process CONSENSUS {
    publishDir "${out}/${sample}", mode: 'copy'
    debug true

    input:
    path(genelist)
    path(db)
    tuple val(sample), path(bam)
    path fasta
    val out

    output:
    path("*.fa"), emit: consensus

    script:
    """
    mkdir ${sample}
    
    bcftools mpileup -Ou -f ${fasta} ${bam} > ${sample}.bcf
    bcftools call --ploidy 1 -m ${sample}.bcf -Oz -o ${sample}.vcf.gz && \\
    tabix ${sample}.vcf.gz && \\

    for gene in \$(cat ${genelist}); do
        COORDS=\$(awk -F'\t' '{print \$5}' <(grep -w "\$gene" ${db}))
        samtools faidx ${fasta} \$COORDS | \\
        bcftools consensus ${sample}.vcf.gz > \$gene.fa
    done
 # for gene in \$(cat ${genelist}); do
 #        COORDS=\$(awk -F'\t' '{print \$4}' <(grep -w "\$gene" ${db}))
 #        
 #        samtools view -b ${bam} \$COORDS | \\
 #        samtools sort -o \$gene.bam && \\
 #        samtools index \$gene.bam && \\
 #        bcftools mpileup -Ou -f ${fasta} \$gene.bam > \$gene.bcf
 #
 #        bcftools call --ploidy 0 -m \$gene.bcf -Oz -o \$gene.vcf.gz && \\
 #        tabix \$gene.vcf.gz && \\
 #        samtools faidx ${fasta} \$COORDS | \\
 #        bcftools consensus \$gene.vcf.gz > \$gene.fa
 #
 #        rm \$gene.bam \$gene.bam.bai \$gene.bcf \$gene.vcf.gz \$gene.vcf.gz.tbi
 #    done
    """
}

workflow {
    fasta = file(params.fasta)
    fastqs = Channel.fromFilePairs("${params.dir}/*_{1,2}.fastq.gz")
    genelist = file(params.genelist)
    db = file(params.db)
    out = file(params.out)
    out.mkdir()

    idx = INDEX(fasta)
    trim = TRIMMOMATIC(fastqs, out)
    bam = DEDUP(idx, trim)
    CONSENSUS(genelist, db, bam, fasta, out)
}
