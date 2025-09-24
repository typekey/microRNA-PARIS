process EXTRACT_UMI {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bin = params.bin
    def new_meta = meta + [single: true]

    """
    $bin/extract_umi.py \\
        --read1 ${reads[0]} \\
        --read2 ${reads[1]} \\
        --output_file ${prefix}.umi.fastq \\

    gzip ${prefix}.umi.fastq
    """
}
