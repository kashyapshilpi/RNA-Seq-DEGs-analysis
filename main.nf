#!/usr/bin/env nextflow
/*
 * RNA-seq differential expression pipeline
 *
 *   FastQC -> Trimmomatic -> FastQC -> HISAT2 -> HTSeq-count -> count matrix
 *          -> DESeq2 + edgeR -> MultiQC
 *
 * Usage:
 *   nextflow run kashyapshilpi/RNA-Seq-DEGs-analysis -profile test,docker
 *   nextflow run kashyapshilpi/RNA-Seq-DEGs-analysis -profile conda \
 *       --input samplesheet.csv --fasta genome.fa --gtf genes.gtf
 */

include { FASTQC as FASTQC_RAW     } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { TRIMMOMATIC              } from './modules/trimmomatic'
include { HISAT2_BUILD             } from './modules/hisat2_build'
include { HISAT2_ALIGN             } from './modules/hisat2_align'
include { HTSEQ_COUNT              } from './modules/htseq_count'
include { MERGE_COUNTS             } from './modules/merge_counts'
include { DESEQ2                   } from './modules/deseq2'
include { EDGER                    } from './modules/edger'
include { MULTIQC                  } from './modules/multiqc'

workflow {

    if (!params.input) { error "Missing --input (samplesheet CSV with columns: sample,fastq_1,fastq_2,condition)" }
    if (!params.fasta) { error "Missing --fasta (reference genome FASTA)" }
    if (!params.gtf)   { error "Missing --gtf (gene annotation GTF)" }

    log.info """
    R N A - S E Q   D E G   P I P E L I N E
    =======================================
    input      : ${params.input}
    fasta      : ${params.fasta}
    gtf        : ${params.gtf}
    contrast   : ${params.contrast_target} vs ${params.contrast_reference}
    outdir     : ${params.outdir}
    """.stripIndent()

    samplesheet = file(params.input, checkIfExists: true)

    ch_reads = channel.fromPath(samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2 || !row.condition) {
                error "Samplesheet rows need sample, fastq_1, fastq_2 and condition: ${row}"
            }
            def meta = [id: row.sample, condition: row.condition]
            [meta, [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]]
        }

    fasta    = file(params.fasta, checkIfExists: true)
    gtf      = file(params.gtf, checkIfExists: true)
    adapters = file(params.adapters, checkIfExists: true)

    // QC and trimming: bad reads cause false alignments and counts
    FASTQC_RAW(ch_reads)
    TRIMMOMATIC(ch_reads, adapters)
    FASTQC_TRIMMED(TRIMMOMATIC.out.reads)

    // Alignment and counting: reads -> genome position -> reads per gene
    HISAT2_BUILD(fasta)
    HISAT2_ALIGN(TRIMMOMATIC.out.reads, HISAT2_BUILD.out.index)
    HTSEQ_COUNT(HISAT2_ALIGN.out.bam, gtf)

    MERGE_COUNTS(HTSEQ_COUNT.out.counts.map { meta, counts -> counts }.collect())

    // Differential expression with two methods: genes found by both are the most reliable
    DESEQ2(MERGE_COUNTS.out.matrix, samplesheet)
    EDGER(MERGE_COUNTS.out.matrix, samplesheet)

    // Combined QC report
    ch_multiqc = FASTQC_RAW.out.reports
        .mix(FASTQC_TRIMMED.out.reports)
        .mix(TRIMMOMATIC.out.log)
        .mix(HISAT2_ALIGN.out.log)
        .mix(HTSEQ_COUNT.out.counts.map { meta, counts -> counts })
        .collect()
    MULTIQC(ch_multiqc)
}
