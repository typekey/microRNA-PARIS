process GENE_ANNOTATION_PEAK_CHIPSEEKER {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.anno.chipseeker.tsv"), emit: annotate_peaks

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"

    """
    run_annotation_chipseeker.R \
        --input_bed_files  ${bed} \
        --output_file ${prefix}.anno.chipseeker.tsv \
        --organism hg38
    """
}
