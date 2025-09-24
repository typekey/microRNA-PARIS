process EXTRACT_POOL_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.mirna_target.fastq.gz"), emit: fastq

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def bin      = params.bin
    def new_meta = meta + [single: true]
    def miRNA    = meta.miRNA
    def mismatch    = params.extract_mismatch

    """
    zcat ${reads[0]} | seqkit grep -j 36 -p ${meta.miRNA} -m ${mismatch} -s >> ${prefix}.${meta.miRNA_name}.${meta.miRNA}.${mismatch}.fastq
    ${bin}/search_paired.py \\
        --input_fastq ${prefix}.${meta.miRNA_name}.${meta.miRNA}.${mismatch}.fastq \\
        --output_alignment ${prefix}.${meta.miRNA_name}.${meta.miRNA}.${mismatch}.mirna_target \\
        --filtered_fastq ${prefix}.${meta.miRNA_name}.${meta.miRNA}.${mismatch}.mirna_target.fastq \\
        --seed_sequence ${miRNA} \\
        --max_mismatches 1 \\
        --prefix_length 7 
    gzip ${prefix}.${meta.miRNA_name}.${meta.miRNA}.${mismatch}.fastq
    gzip ${prefix}.${meta.miRNA_name}.${meta.miRNA}.${mismatch}.mirna_target
    gzip ${prefix}.${meta.miRNA_name}.${meta.miRNA}.${mismatch}.mirna_target.fastq
    """
}
