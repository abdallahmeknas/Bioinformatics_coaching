#!/usr/bin/env nextflow
/*
 * Mini variant calling pipeline for teaching
 *
 *   raw reads -> FastQC -> fastp -> BWA -> samtools -> bcftools -> MultiQC
 *
 * Run it:            nextflow run main.nf
 * Run it again:      nextflow run main.nf -resume
 * With conda:        nextflow run main.nf -profile conda
 * With containers:   nextflow run main.nf -profile docker
 */

// ---------------------------------------------------------------------------
// Processes: each one is a single step with its own inputs and outputs
// ---------------------------------------------------------------------------

process FASTQC {
    tag "${sample}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${r1} ${r2}
    """
}

process FASTP {
    tag "${sample}"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.{html,json}'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_R1.trim.fastq.gz"), path("${sample}_R2.trim.fastq.gz"), emit: reads
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.html", emit: html

    script:
    """
    fastp \\
        --in1 ${r1} --in2 ${r2} \\
        --out1 ${sample}_R1.trim.fastq.gz --out2 ${sample}_R2.trim.fastq.gz \\
        --detect_adapter_for_pe \\
        --json ${sample}.fastp.json --html ${sample}.fastp.html \\
        --thread ${task.cpus}
    """
}

process BWA_INDEX {
    tag "${fasta.simpleName}"

    input:
    path fasta

    output:
    path "bwa_index", emit: index

    script:
    """
    mkdir bwa_index
    cp ${fasta} bwa_index/genome.fa
    bwa index bwa_index/genome.fa
    """
}

process BWA_MEM {
    tag "${sample}"

    input:
    tuple val(sample), path(r1), path(r2)
    path index

    output:
    tuple val(sample), path("${sample}.sam"), emit: sam

    script:
    """
    bwa mem \\
        -t ${task.cpus} \\
        -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' \\
        ${index}/genome.fa ${r1} ${r2} > ${sample}.sam
    """
}

process SAMTOOLS_SORT {
    tag "${sample}"
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam ${sam}
    samtools index ${sample}.sorted.bam
    """
}

process SAMTOOLS_STATS {
    tag "${sample}"
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path "${sample}.{flagstat,stats}", emit: reports

    script:
    """
    samtools flagstat ${bam} > ${sample}.flagstat
    samtools stats ${bam} > ${sample}.stats
    """
}

process SAMTOOLS_FAIDX {
    tag "${fasta.simpleName}"

    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai

    script:
    """
    samtools faidx ${fasta}
    """
}

process BCFTOOLS_CALL {
    tag "cohort"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path bams
    path bais
    path fasta
    path fai

    output:
    path "cohort.vcf.gz", emit: vcf
    path "cohort.vcf.gz.csi", emit: index
    path "cohort.bcftools_stats.txt", emit: stats

    script:
    """
    bcftools mpileup --threads ${task.cpus} -f ${fasta} -a FORMAT/AD,FORMAT/DP -Ou ${bams} \\
        | bcftools call -mv -Oz -o cohort.vcf.gz
    bcftools index cohort.vcf.gz
    bcftools stats cohort.vcf.gz > cohort.bcftools_stats.txt
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html", emit: report
    path "*_data", emit: data

    script:
    """
    multiqc --force --filename multiqc_report.html .
    """
}

// ---------------------------------------------------------------------------
// Workflow: wires the processes together with channels
// ---------------------------------------------------------------------------

workflow {
    log.info """
    ==============================================
     Mini variant calling pipeline
     samplesheet : ${params.samplesheet}
     reference   : ${params.reference}
     outdir      : ${params.outdir}
    ==============================================
    """.stripIndent()

    // One item per sample: [ sample, read1, read2 ]
    reads_ch = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)) }

    reference = file(params.reference, checkIfExists: true)

    // QC and trimming
    FASTQC(reads_ch)
    FASTP(reads_ch)

    // Alignment
    BWA_INDEX(reference)
    BWA_MEM(FASTP.out.reads, BWA_INDEX.out.index)
    SAMTOOLS_SORT(BWA_MEM.out.sam)
    SAMTOOLS_STATS(SAMTOOLS_SORT.out.bam)

    // Joint variant calling across all samples
    SAMTOOLS_FAIDX(reference)
    bams = SAMTOOLS_SORT.out.bam.map { sample, bam, bai -> bam }.collect(sort: true)
    bais = SAMTOOLS_SORT.out.bam.map { sample, bam, bai -> bai }.collect(sort: true)
    BCFTOOLS_CALL(bams, bais, reference, SAMTOOLS_FAIDX.out.fai)

    // One report to rule them all
    qc_files = FASTQC.out.reports
        .mix(FASTP.out.json, SAMTOOLS_STATS.out.reports, BCFTOOLS_CALL.out.stats)
        .collect()
    MULTIQC(qc_files)
}
