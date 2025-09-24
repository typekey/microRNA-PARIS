#!/usr/bin/env nextflow
/*
----------------------------------------------------------------------------------------
File    :   main.nf
Time    :   2025/09/23 17:36:31
Author  :   Lei Zheng
Version :   1.0
Contact :   type.zheng@gmail.com
Github  :   https://github.com/typekey
----------------------------------------------------------------------------------------
*/


include { MACS3_CALLPEAK } from '../../../modules/local/macs3'

workflow CALL_PEAK {
    take:
    ch_genome_bam       // channel: [ val(meta), path(bam) ]
    with_control        // boolean: whether to use control samples

    main:
    if (params.with_control) {
        ch_genome_bam
            .map {
                meta, bam ->
                    meta.control ? null : [ meta.id, bam ]
            }
            .set { ch_bam_merged_control }
        
        ch_genome_bam
            .map {
                meta, bam ->
                    control_id = meta.control + "_" + meta.control_replicate
                    meta.control ? [ control_id, meta, bam ] : null
            }
            .combine( ch_bam_merged_control, by: 0 )

            .map { it -> [ it[1] , it[2], it[3] ] }
            .set { ch_bam_replicate }
    } else {
        ch_genome_bam
            .map {
                meta, bam ->
                    [ meta , bam, [] ]
            }
            .set { ch_bam_replicate }
    }

    MACS3_CALLPEAK (
        ch_bam_replicate
    )
    // -------------------
    emit:
    bed = MACS3_CALLPEAK.out.bed
}
