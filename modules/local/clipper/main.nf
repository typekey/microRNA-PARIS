process CLIPPER {
    tag "$meta.id"
    label 'process_medium'

    container 'docker://brianyee/clipper:6594e71'
    // container 'docker://brianyee/clipper:5d865bb'

    input:
    tuple val(meta), path(bam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo "Running CLIPPER"
    clipper --bam $bam --species GRCh38 --outfile ${prefix}.peaks.bed
    """
}
