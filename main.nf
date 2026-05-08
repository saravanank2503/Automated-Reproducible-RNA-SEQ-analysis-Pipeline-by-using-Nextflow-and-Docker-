#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.samplesheet   = "$projectDir/samplesheet.csv"
params.outdir        = "$projectDir/results"
params.star_index    = null
params.gtf           = null
params.condition_col = "condition"
params.deseq2_fdr    = 0.05
params.deseq2_lfc    = 1.0
params.max_memory    = "48.GB"
params.max_cpus      = 6
params.max_time      = "24.h"

log.info """
==============================================
  RNA-SEQ TEST PIPELINE
  Samples  : ${params.samplesheet}
  STAR idx : ${params.star_index}
  GTF      : ${params.gtf}
  Output   : ${params.outdir}
==============================================
""".stripIndent()

if (!params.star_index) error "ERROR: --star_index is required"
if (!params.gtf)        error "ERROR: --gtf is required"

include { FASTQC as FASTQC_RAW  } from './modules/fastqc'
include { TRIMGALORE             } from './modules/trimgalore'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc'
include { STAR_ALIGN             } from './modules/star_align'
include { SAMTOOLS_SORT_INDEX    } from './modules/samtools'
include { FEATURECOUNTS          } from './modules/featurecounts'
include { MERGE_COUNTS           } from './modules/merge_counts'
include { DESEQ2_ANALYSIS        } from './modules/deseq2'
include { MULTIQC                } from './modules/multiqc'

def parse_samplesheet(csv) {
    Channel.fromPath(csv)
        .splitCsv(header: true, strip: true)
        .map { row ->
            ['sample','fastq_1','condition'].each { col ->
                if (!row[col]) error "Samplesheet missing: '${col}' in ${row}"
            }
            def meta = [
                id        : row.sample,
                condition : row.condition,
                replicate : row.replicate ?: '1',
                batch     : row.batch     ?: 'batch1',
                paired    : row.fastq_2   ? true : false
            ]
            def reads = meta.paired
                ? [ file(row.fastq_1, checkIfExists: true),
                    file(row.fastq_2, checkIfExists: true) ]
                : [ file(row.fastq_1, checkIfExists: true) ]
            return [ meta, reads ]
        }
}

workflow {

    ch_reads      = parse_samplesheet(params.samplesheet)
    ch_star_index = Channel.value(file(params.star_index))
    ch_gtf        = Channel.value(file(params.gtf))

    FASTQC_RAW   ( ch_reads )
    TRIMGALORE   ( ch_reads )
    FASTQC_TRIM  ( TRIMGALORE.out.reads )

    STAR_ALIGN        ( TRIMGALORE.out.reads, ch_star_index, ch_gtf )
    SAMTOOLS_SORT_INDEX ( STAR_ALIGN.out.bam )

    FEATURECOUNTS (
        SAMTOOLS_SORT_INDEX.out.bam,
        ch_gtf
    )

    MERGE_COUNTS ( FEATURECOUNTS.out.counts.map { meta, counts -> counts }.collect() )

    DESEQ2_ANALYSIS (
        MERGE_COUNTS.out.matrix,
        file(params.samplesheet),
        params.condition_col,
        params.deseq2_fdr,
        params.deseq2_lfc
    )

    // ── Strip meta maps before collecting for MultiQC ─────────────────────
    ch_qc = Channel.empty()
        .mix( FASTQC_RAW.out.zip.map   { meta, zip -> zip }.flatten()  )
        .mix( FASTQC_TRIM.out.zip.map  { meta, zip -> zip }.flatten()  )
        .mix( TRIMGALORE.out.log.flatten()                              )
        .mix( STAR_ALIGN.out.log.flatten()                              )
        .mix( FEATURECOUNTS.out.summary.flatten()                       )
        .collect()

    MULTIQC ( ch_qc )
}

workflow.onComplete {
    log.info workflow.success
        ? "PIPELINE COMPLETE → ${params.outdir}"
        : "PIPELINE FAILED  → nextflow log ${workflow.runName} -f name,status,exit,workdir"
}
