//
// Alignment with STAR
//

include { MACS3_CALLPEAK as CALL_PEAK_PLUS } from '../../../modules/local/macs3_strand'
include { MACS3_CALLPEAK as CALL_PEAK_MINUS } from '../../../modules/local/macs3_strand'
include { SPLIT_STRAND } from '../../../modules/local/split_strand'
include { MERGE_STRAND } from '../../../modules/local/merge_strand'

workflow CALL_PEAK_STRAND {
    take:
    ch_genome_bam       // channel: [ val(meta), path(bam) ]
    with_control        // boolean: whether to use control samples

    main:
    SPLIT_STRAND (
        ch_genome_bam
    )

    SPLIT_STRAND.out.plus_bam
        .map {
            meta, bam ->
                [ meta, bam ]
        }
        .set { ch_genome_bam_plus }

    if (params.with_control) {
        ch_genome_bam_plus
            .map {
                meta, bam ->
                    meta.control ? null : [ meta.id, bam ]
            }
            .set { ch_bam_merged_control }
        
        ch_genome_bam_plus
            .map {
                meta, bam ->
                    control_id = meta.control + "_" + meta.control_replicate
                    meta.control ? [ control_id, meta, bam ] : null
            }
            .combine( ch_bam_merged_control, by: 0 )

            .map { it -> [ it[1] , it[2], it[3] ] }
            .set { ch_bam_replicate_plus }
    } else {
        ch_genome_bam_plus
            .map {
                meta, bam ->
                    [ meta , bam, [] ]
            }
            .set { ch_bam_replicate_plus }
    }


    SPLIT_STRAND.out.minus_bam
        .map {
            meta, bam ->
                [ meta, bam ]
        }
        .set { ch_genome_bam_minus }

    if (params.with_control) {
        ch_genome_bam_minus
            .map {
                meta, bam ->
                    meta.control ? null : [ meta.id, bam ]
            }
            .set { ch_bam_merged_control }
        
        ch_genome_bam_minus
            .map {
                meta, bam ->
                    control_id = meta.control + "_" + meta.control_replicate
                    meta.control ? [ control_id, meta, bam ] : null
            }
            .combine( ch_bam_merged_control, by: 0 )

            .map { it -> [ it[1] , it[2], it[3] ] }
            .set { ch_bam_replicate_minus }
    } else {
        ch_genome_bam_minus
            .map {
                meta, bam ->
                    [ meta , bam, [] ]
            }
            .set { ch_bam_replicate_minus }
    }
    
    // ch_bam_replicate_minus.view()
    // ch_genome_bam_minus.view()

    CALL_PEAK_PLUS (
        ch_bam_replicate_plus,
        "plus"
    )

    CALL_PEAK_MINUS (
        ch_bam_replicate_minus,
        "minus"
    )

    CALL_PEAK_PLUS.out.peak
        .combine( CALL_PEAK_MINUS.out.peak, by: 0 )
        .set { ch_peak }

    MERGE_STRAND (
        ch_peak
    )
    
    emit:
    bed = MERGE_STRAND.out.bed
}
