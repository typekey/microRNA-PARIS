process GENE_ANNOTATION_PEAK {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.annotatePeaks.txt"), emit: annotate_peaks

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"

    """
    run_annotation_peak.py -s hg38 --bed_file ${bed} --threads $task.cpus --output ${prefix}.annotatePeaks.txt
    """
}
