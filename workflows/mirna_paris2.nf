#!/usr/bin/env nextflow
/*
----------------------------------------------------------------------------------------
File    :   mirna_paris2.nf
Time    :   2025/01/07 22:50:26
Author  :   Lei Zheng
Version :   1.0
Contact :   type.zheng@gmail.com
Github  :   https://github.com/typekey
----------------------------------------------------------------------------------------
*/

include { CAT_FASTQ                      } from '../modules/local/cat/fastq'
include { FROM_CSV                       } from '../modules/local/import_data/from_csv'
include { READ_MIRNA_TABLE               } from '../modules/local/import_data/read_mirna_table'
include { READ_POOL_TABLE                } from '../modules/local/import_data/read_pool_table'
include { MIRNA_PAIRED                   } from '../modules/local/mirna_paired'
include { BEDTOOLS_CALL_PEAK             } from '../modules/local/bedtools_call_peak'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions'

include { FASTQ_FASTQC_EXTRACT_UMI_TRIMGALORE_PARIED2SINGLE } from '../subworkflows/local/fastq_fastqc_extract_umi_trimgalore_paired2single'
include { FASTQ_FASTQC_UMITOOLS_CUTADAPT_paired2single } from '../subworkflows/local/fastq_fastqc_extract_umi_cutadapt_paired2single'
include { FASTQ_FASTQC_UMITOOLS_CUTADAPT } from '../subworkflows/nf-core/fastq_fastqc_umitools_cutadapt'
include { PLOT_PEAK_NUMBER            } from '../modules/local/plot/peak_number'
include { GENE_ANNOTATION_PEAK            } from '../modules/local/annotation/gene_annotation_peak'
include { GENE_ANNOTATION_PEAK_CHIPSEEKER            } from '../modules/local/annotation/gene_annotation_peak_chipseeker'
include { PLOT_PEAK_ANNOTATION_CHIPSEEKER            } from '../modules/local/plot/peak_annotation_chipseeker'

include { BBMAP_BBSPLIT               } from '../modules/nf-core/bbmap/bbsplit'
include { PREPARE_GENOME                                    } from '../subworkflows/local/prepare_genome'
include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../subworkflows/nf-core/fastq_subsample_fq_salmon'
include { ALIGN_STAR as ALIGN_REPEAT_START                                        } from '../subworkflows/local/align_star'
include { ALIGN_STAR as ALIGN_GENE_START                                          } from '../subworkflows/local/align_star'
include { FEATRURECOUNTS_GENE                                                  } from '../modules/local/featurecount_gene'
include { FEATRURECOUNTS_REPEAT                                                   } from '../modules/local/featurecount_repeat'
include { EXTRACT_FASTQ                                                           } from '../modules/local/extract_fastq'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENE        } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_REPEAT        } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BAM_MARKDUPLICATES_PICARD        } from '../subworkflows/nf-core/bam_markduplicates_picard'
include { MULTIQC_CUSTOM_PEAKS               } from '../modules/local/multiqc/custom_modules/frip_score_multiqc'

include { SAMTOOLS_SORT               } from '../modules/nf-core/samtools/sort'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../modules/local/umitools_prepareforrsem'
include { BAM_SORT_STATS_SAMTOOLS          } from '../subworkflows/nf-core/bam_sort_stats_samtools'
include { CLIPPER               } from '../modules/local/clipper'
include { PURECLIP               } from '../modules/local/pureclip'
include { SAMTOOLS_FILTER_CHROM_NORM        } from '../modules/nf-core/samtools/filter_norm'
include { CALL_PEAK as RAW_CALL_PEAK                         } from '../subworkflows/local/call_peak'
include { CALL_PEAK as DEDUP_CALL_PEAK                        } from '../subworkflows/local/call_peak'
include { CALL_PEAK_STRAND as RAW_CALL_PEAK_STRAND                        } from '../subworkflows/local/call_peak_strand'
include { CALL_PEAK_STRAND as DEDUP_CALL_PEAK_STRAND                        } from '../subworkflows/local/call_peak_strand'

include { FRIP_SCORE               } from '../modules/local/frip_score'

include { MULTIQC                            } from '../modules/local/multiqc'
include { HOMER_MOTIF                                       } from '../modules/local/homer/motif'
include { HOMER_ANNOTATEPEAKS                               } from '../modules/local/homer/annotatepeaks'
include { FASTQSCREEN                                         } from '../modules/local/fastq_screen'

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def summary_params = paramsSummaryMap(workflow)

// --------------------------------------------------------------------
// Print parameter summary log to screen
// --------------------------------------------------------------------
log.info paramsSummaryLog(workflow)

workflow miRNA_PARIS2 {
    println("Starting PARIS2 workflow...")
    ch_versions = Channel.empty()

    // -----------------------------------------------------------------------------------------------------------
    // Uncompress and prepare reference genome files
    // -----------------------------------------------------------------------------------------------------------
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    def prepareToolIndices  = []
    if (!params.skip_alignment) { prepareToolIndices << params.aligner }

    def filterGtf =
        ((
            !params.skip_alignment && params.aligner
        ) ||
        (
            !params.skip_pseudo_alignment && params.pseudo_aligner
        ) ||
        (
            !params.transcript_fasta
        )) &&
        (
            !params.skip_gtf_filter
        )

    def is_aws_igenome = false
    if (params.fasta && params.gtf) {
        if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
            is_aws_igenome = true
        }
    }

    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.gene_bed,
        params.splicesites,
        params.star_index,
        params.hisat2_index,
        params.gencode,
        is_aws_igenome,
        biotype,
        prepareToolIndices,
        filterGtf
    )


    // -----------------------------------------------------------------------------------------------------------
    // Create input channel from input file
    // -----------------------------------------------------------------------------------------------------------

    fq_ch = FROM_CSV(params.input)

    CAT_FASTQ (
        fq_ch
    )
    .set { ch_cat_fastq }

    ch_filtered_reads      = Channel.empty()
    cutadapt_json          = Channel.empty()
    umi_extract_log        = Channel.empty()
    ch_fastqc_raw_multiqc  = Channel.empty()
    ch_fastqc_trim_multiqc = Channel.empty()
    ch_trim_log_multiqc    = Channel.empty()
    ch_trim_read_count     = Channel.empty()

    if (params.trimmer == 'cutadapt') {
        FASTQ_FASTQC_UMITOOLS_CUTADAPT_paired2single (
            fq_ch,
            params.skip_fastqc,
            params.with_umi,
            params.skip_umi_extract,
            params.skip_trimming,
            params.umi_discard_read,
            params.min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_CUTADAPT_paired2single.out.reads
        cutadapt_json      = FASTQ_FASTQC_UMITOOLS_CUTADAPT_paired2single.out.cutadapt_json
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_CUTADAPT_paired2single.out.fastqc_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_CUTADAPT_paired2single.out.trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_CUTADAPT_paired2single.out.trim_log
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_CUTADAPT_paired2single.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_CUTADAPT_paired2single.out.versions)
    }

    ch_fastq_screen_multiqc = Channel.empty()
    if (!params.skip_fastq_screen) {
        FASTQSCREEN (
            ch_filtered_reads,
            params.fastq_screen_conf
        )
        ch_fastq_screen_multiqc = FASTQSCREEN.out.txt
    }

    MIRNA_PAIRED (
        ch_filtered_reads
    )

    // -------------------------------------------------------------------------
    // Alignment 
    // -------------------------------------------------------------------------

    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    umi_dedup_log           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star') {
        if (params.map_target == 'gene') {
            ALIGN_GENE_START (
                ch_filtered_reads,
                PREPARE_GENOME.out.star_index.map { [ [:], it ] },
                PREPARE_GENOME.out.gtf.map { [ [:], it ] },
                params.star_ignore_sjdbgtf,
                '',
                params.seq_center ?: '',
                is_aws_igenome,
                PREPARE_GENOME.out.fasta.map { [ [:], it ] },
                "gene"
            )
            ch_gene_bam        = ALIGN_GENE_START.out.bam
            ch_gene_bam_index  = ALIGN_GENE_START.out.bai
            ch_gene_unmapped_1 = ALIGN_GENE_START.out.unmapped_1
            ch_gene_transcriptome_bam = ALIGN_GENE_START.out.bam_transcript
            ch_gene_samtools_stats    = ALIGN_GENE_START.out.stats
            ch_gene_samtools_flagstat = ALIGN_GENE_START.out.flagstat
            ch_gene_samtools_idxstats = ALIGN_GENE_START.out.idxstats
            ch_gene_star_multiqc      = ALIGN_GENE_START.out.log_final
            if (params.bam_csi_index) {
                ch_repeat_bam_index = ALIGN_GENE_START.out.csi
            }
            ch_versions = ch_versions.mix(ALIGN_GENE_START.out.versions)

        } else if (params.map_target == 'repeat') {
            ALIGN_REPEAT_START (
                ch_filtered_reads,
                PREPARE_GENOME.out.star_index.map { [ [:], it ] },
                PREPARE_GENOME.out.gtf.map { [ [:], it ] },
                params.star_ignore_sjdbgtf,
                '',
                params.seq_center ?: '',
                is_aws_igenome,
                PREPARE_GENOME.out.fasta.map { [ [:], it ] },
                "repeat"
            )
            ch_repeat_bam        = ALIGN_REPEAT_START.out.bam
            ch_repeat_bam_index  = ALIGN_REPEAT_START.out.bai
            ch_repeat_unmapped_1 = ALIGN_REPEAT_START.out.unmapped_1
            ch_repeat_transcriptome_bam = ALIGN_REPEAT_START.out.bam_transcript
            ch_repeat_samtools_stats    = ALIGN_REPEAT_START.out.stats
            ch_repeat_samtools_flagstat = ALIGN_REPEAT_START.out.flagstat
            ch_repeat_samtools_idxstats = ALIGN_REPEAT_START.out.idxstats
            ch_repeat_star_multiqc      = ALIGN_REPEAT_START.out.log_final
            if (params.bam_csi_index) {
                ch_repeat_bam_index = ALIGN_REPEAT_START.out.csi
            }

        } else if (params.map_target == 'both') {
            ALIGN_REPEAT_START (
                ch_filtered_reads,
                PREPARE_GENOME.out.star_index.map { [ [:], it ] },
                PREPARE_GENOME.out.gtf.map { [ [:], it ] },
                params.star_ignore_sjdbgtf,
                '',
                params.seq_center ?: '',
                is_aws_igenome,
                PREPARE_GENOME.out.fasta.map { [ [:], it ] },
                "repeat"
            )
            ch_repeat_bam        = ALIGN_REPEAT_START.out.bam
            ch_repeat_bam_index  = ALIGN_REPEAT_START.out.bai
            ch_repeat_unmapped_1 = ALIGN_REPEAT_START.out.unmapped_1
            ch_repeat_transcriptome_bam = ALIGN_REPEAT_START.out.bam_transcript
            ch_repeat_samtools_stats    = ALIGN_REPEAT_START.out.stats
            ch_repeat_samtools_flagstat = ALIGN_REPEAT_START.out.flagstat
            ch_repeat_samtools_idxstats = ALIGN_REPEAT_START.out.idxstats
            ch_repeat_star_multiqc      = ALIGN_REPEAT_START.out.log_final
            if (params.bam_csi_index) {
                ch_repeat_bam_index = ALIGN_REPEAT_START.out.csi
            }
            ch_versions = ch_versions.mix(ALIGN_REPEAT_START.out.versions)

            ALIGN_GENE_START (
                ch_repeat_unmapped_1,
                PREPARE_GENOME.out.star_index.map { [ [:], it ] },
                PREPARE_GENOME.out.gtf.map { [ [:], it ] },
                params.star_ignore_sjdbgtf,
                '',
                params.seq_center ?: '',
                is_aws_igenome,
                PREPARE_GENOME.out.fasta.map { [ [:], it ] },
                "gene"
            )
            ch_gene_bam        = ALIGN_GENE_START.out.bam
            ch_gene_bam_index  = ALIGN_GENE_START.out.bai
            ch_gene_unmapped_1 = ALIGN_GENE_START.out.unmapped_1
            ch_gene_transcriptome_bam = ALIGN_GENE_START.out.bam_transcript
            ch_gene_samtools_stats    = ALIGN_GENE_START.out.stats
            ch_gene_samtools_flagstat = ALIGN_GENE_START.out.flagstat
            ch_gene_samtools_idxstats = ALIGN_GENE_START.out.idxstats
            ch_gene_star_multiqc      = ALIGN_GENE_START.out.log_final
            if (params.bam_csi_index) {
                ch_repeat_bam_index = ALIGN_GENE_START.out.csi
            }
            ch_versions = ch_versions.mix(ALIGN_GENE_START.out.versions)

        } else {
            error "Unknown mapping target: ${params.map_target}"
        }
    }

    // -------------------------------------------------------------------
    // Remove duplicate reads from BAM file based on UMIs
    // -------------------------------------------------------------------

    if (params.with_umi) {
        // Deduplicate genome BAM file before downstream analysis
        if (params.map_target == 'gene' || params.map_target == 'both') {
            BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENE (
                ch_gene_bam.join(ch_gene_bam_index, by: [0]),
                params.umitools_dedup_stats
            )
            ch_gene_bam        = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENE.out.bam
            ch_gene_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENE.out.bai
            ch_gene_umi_dedup_log  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENE.out.umi_dedup_log
            ch_gene_samtools_stats    = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENE.out.stats
            ch_gene_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENE.out.flagstat
            ch_gene_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENE.out.idxstats
            if (params.bam_csi_index) {
                ch_gene_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENE.out.csi
            }
            ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENE.out.versions)
            
            // SAMTOOLS_FILTER_CHROM_NORM (
            //     ch_genome_bam
            // )
        } 
        
        if (params.map_target == 'repeat' || params.map_target == 'both') {
            BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_REPEAT (
                ch_repeat_bam.join(ch_repeat_bam_index, by: [0]),
                params.umitools_dedup_stats
            )
            ch_repeat_bam        = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_REPEAT.out.bam
            ch_repeat_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_REPEAT.out.bai
            ch_repeat_umi_dedup_log  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_REPEAT.out.umi_dedup_log
            ch_repeat_samtools_stats    = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_REPEAT.out.stats
            ch_repeat_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_REPEAT.out.flagstat
            ch_repeat_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_REPEAT.out.idxstats
            if (params.bam_csi_index) {
                ch_repeat_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_REPEAT.out.csi
            }
            ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_REPEAT.out.versions)
            
            // SAMTOOLS_FILTER_CHROM_NORM (
            //     ch_genome_bam
            // )
        }
    }

    // -------------------------------------------------------------------
    // Call peaks and annotate peaks
    // -------------------------------------------------------------------

    DEDUP_CALL_PEAK_STRAND (
        ch_gene_bam,
        params.with_control
    )

    GENE_ANNOTATION_PEAK (
        DEDUP_CALL_PEAK_STRAND.out.bed,
    )

    GENE_ANNOTATION_PEAK_CHIPSEEKER (
        DEDUP_CALL_PEAK_STRAND.out.bed,
    )

    // BEDTOOLS_CALL_PEAK(ch_gene_bam)

    // GENE_ANNOTATION_PEAK (
    //     BEDTOOLS_CALL_PEAK.out.bed,
    // )

    // GENE_ANNOTATION_PEAK_CHIPSEEKER (
    //     BEDTOOLS_CALL_PEAK.out.bed,
    // )
}
