process FEATRURECOUNTS_GENE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.count"), emit: count
    tuple val(meta), path("*annotation.txt"), emit: annotation

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def bin      = params.bin
    def new_meta = meta + [single: true]
    def miRNA    = meta.miRNA
    def miRNA_name    = meta.miRNA_name
    def mismatch    = params.extract_mismatch

    """
    run_featurecounts.py -s hg38  --gene_type gene -b ${bam} --output ${miRNA_name}_featurecount_result_gene.count
    ${bin}/annotation.py --count_file ${miRNA_name}_featurecount_result_gene.count --miRNA ${miRNA_name} --output ${miRNA_name}_gene_annotation.txt

    """
}
