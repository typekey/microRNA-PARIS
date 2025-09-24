//
// Read QC, UMI extraction and trimming
//

include { FASTQC                         } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as TRIM_FASTQC          } from '../../../modules/nf-core/fastqc/main'
include { CUTADAPT                       } from '../../../modules/local/cutadapt/main'
include { EXTRACT_UMI                    } from '../../../modules/local/extract_umi'
include { READ_MIRNA_TABLE               } from '../../../modules/local/import_data/read_mirna_table'
include { READ_POOL_TABLE                } from '../../../modules/local/import_data/read_pool_table'
include { EXTRACT_POOL_FASTQ             } from '../../../modules/local/extract_pool_fastq'

workflow FASTQ_FASTQC_UMITOOLS_CUTADAPT_paired2single {
    take:
    reads             // channel: [ val(meta), [ reads ] ]
    skip_fastqc       // boolean: true/false
    with_umi          // boolean: true/false
    skip_umi_extract  // boolean: true/false
    skip_trimming     // boolean: true/false
    umi_discard_read  // integer: 0, 1 or 2
    min_trimmed_reads // integer: > 0

    main:
    ch_versions = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()
    if (!skip_fastqc) {
        FASTQC (reads)
        fastqc_html = FASTQC.out.html
        fastqc_zip  = FASTQC.out.zip
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    umi_reads = reads
    umi_extract_log   = Channel.empty()
    if (!skip_umi_extract) {
        EXTRACT_UMI (
            umi_reads
        )
        .set { umi_reads }
    }
    // umi_reads.view()

    if (params.is_pooled) {
        READ_MIRNA_TABLE(
            params.mirna_table
        )

        // umi_reads.view()
        // umi_reads
        //     .combine(READ_MIRNA_TABLE.out.mirna_ch)
        //     .map { meta, reads, mirna_info ->
        //         def new_meta = meta + [
        //             mirna_info: mirna_info,
        //             miRNA     : mirna_info.extract_seq,
        //             miRNA_name: mirna_info.miRNA,
        //             seed      : mirna_info.seed,
        //             id        : "${meta.id}_${mirna_info.miRNA}",
        //         ]
        //         return [new_meta, reads]
        //     }
        //     .set { combined_mirna_reads_ch }

         combined_mirna_reads_ch = Channel.empty()

        if (params.mirna_table) {
            READ_MIRNA_TABLE(
                params.mirna_table
            )
            // CAT_FASTQ.out.reads.view()
            // READ_MIRNA_TABLE.out.mirna_ch.view()

            umi_reads
                .map { meta, reads ->
                    tuple( meta.miRNA, meta, reads )
                }
                .set { keyed_reads_ch }

            // keyed_reads_ch.view()

            READ_MIRNA_TABLE.out.mirna_ch
                .map { meta ->
                    // 把 meta.miRNA 当作第一个元素做 key
                    tuple( meta.miRNA, meta )
                }
                .set { keyed_mirna_ch }

            // keyed_mirna_ch.view()

            keyed_reads_ch
                .join( keyed_mirna_ch, by: [0,0] )
                .map { miRNA_key, meta, reads, mirna_meta ->
                    def new_meta = meta + [
                    // mirna_info:  mirna_meta,
                        miRNA        : mirna_meta.extract_seq,
                        miRNA_name   : mirna_meta.miRNA,
                        seed         : mirna_meta.seed,
                        target       : mirna_meta.target,
                        extract_seq  : mirna_meta.extract_seq,
                        single_end   : true,
                        id           : "${meta.id}_${mirna_meta.miRNA}"
                    ]
                    tuple( new_meta, reads )
                }
                .set { combined_single_mirna_reads_ch }
            combined_mirna_reads_ch = combined_mirna_reads_ch.mix(combined_single_mirna_reads_ch)
        }

        if (params.pool_table) {
            READ_POOL_TABLE(
                params.pool_table
            )

            // ch_cat_fastq.view()
            umi_reads
                .combine(READ_POOL_TABLE.out.mirna_ch)
                .map { meta, reads, mirna_info ->
                    if (meta.miRNA==null) {
                        def new_meta = meta + [
                            // mirna_info: mirna_info,
                            miRNA        : mirna_meta.extract_seq,
                            miRNA_name   : mirna_meta.miRNA,
                            seed         : mirna_info.seed,
                            target       : mirna_info.target,
                            extract_seq  : mirna_info.extract_seq,
                            single_end   : true,
                            id           : "${meta.id}_${mirna_info.miRNA}"
                        ]
                        return [new_meta, reads]
                    } 
                }
                .set { combined_pool_mirna_reads_ch }
            combined_mirna_reads_ch = combined_mirna_reads_ch.mix(combined_pool_mirna_reads_ch)
        }

        combined_mirna_reads_ch.view()

        EXTRACT_POOL_FASTQ (
            combined_mirna_reads_ch
        )
        // umi_reads = EXTRACT_POOL_FASTQ.out.fastq
        umi_reads = EXTRACT_POOL_FASTQ.out.fastq
                    .filter { meta, fq -> fq.size() > 10000 }

    } else {
        EXTRACT_FASTQ (
            fq_ch,
        )
        umi_reads = EXTRACT_FASTQ.out.fastq
    }

    trim_reads      = umi_reads
    trim_unpaired   = Channel.empty()
    trim_html       = Channel.empty()
    trim_zip        = Channel.empty()
    trim_log        = Channel.empty()
    trim_read_count = Channel.empty()
    if (!skip_trimming) {
        // cutadapt
        CUTADAPT (umi_reads)
        trim_reads = CUTADAPT.out.reads
        cutadapt_json = CUTADAPT.out.json

        TRIM_FASTQC (trim_reads)
        trim_html = TRIM_FASTQC.out.html
        trim_zip  = TRIM_FASTQC.out.zip
        ch_versions   = ch_versions.mix(CUTADAPT.out.versions.first())
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]
    cutadapt_json

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]

    umi_extract_log            // channel: [ val(meta), [ log ] ]

    trim_unpaired      // channel: [ val(meta), [ reads ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]
    trim_read_count    // channel: [ val(meta), val(count) ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
