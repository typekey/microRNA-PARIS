process HOMER_ANNOTATEPEAKS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bed)
    val(species)

    output:
    tuple val(meta), path("*annotatePeaks.txt"), emit: txt
    tuple val(meta), path("*annStats.txt"), emit: stats, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_annotation_bedtools_carna_by_peak.py -s $species --input_bed ${bed} --output_path carna_anno.tsv
    """
}