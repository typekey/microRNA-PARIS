process MIRNA_PAIRED {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.alignment.gz"), emit: mirna_paired
    tuple val(meta), path("*.motif.gz"), emit: mirna_paired_motif

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def bin      = params.bin
    def new_meta = meta + [single: true]
    def miRNA    = meta.miRNA
    def mismatch    = params.extract_mismatch

    """
    ${bin}/build_mirna_paired.py --input_fastq ${reads[0]} --output_dir ./ --seed_sequence ${meta.seed}
    """
}
