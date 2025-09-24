process FEATRURECOUNTS_REPEAT {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.count"), emit: count
    tuple val(meta), path("*annotation.txt"), emit: annotation

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def bin        = params.bin
    def new_meta   = meta + [single: true]
    def miRNA      = meta.miRNA
    def miRNA_name = meta.miRNA_name
    def mismatch   = params.extract_mismatch

    """
    run_featurecounts.py -s hg38  --gene_type repeat_unique -b ${bam} --output ${miRNA_name}_featurecount_result_repeat.count
    ${bin}/annotation_carna.py --count_file ${miRNA_name}_featurecount_result_repeat.count --miRNA ${miRNA_name} --biotype 'repeat' --output ${miRNA_name}_repeat_annotation.txt

    run_featurecounts.py -s hg38  --gene_type enhancer -b ${bam} --output ${miRNA_name}_featurecount_result_enhancer.count
    ${bin}/annotation_carna.py --count_file ${miRNA_name}_featurecount_result_enhancer.count --miRNA ${miRNA_name} --biotype 'enhancer' --output ${miRNA_name}_enhancer_annotation.txt

    run_featurecounts.py -s hg38  --gene_type promoter_minus -b ${bam} --output ${miRNA_name}_featurecount_result_promoter_minus.count
    ${bin}/annotation_carna.py --count_file ${miRNA_name}_featurecount_result_promoter_minus.count --miRNA ${miRNA_name} --biotype 'promoter_minus' --output ${miRNA_name}_promoter_minus_annotation.txt

    run_featurecounts.py -s hg38  --gene_type promoter_plus -b ${bam} --output ${miRNA_name}_featurecount_result_promoter_plus.count
    ${bin}/annotation_carna.py --count_file ${miRNA_name}_featurecount_result_promoter_plus.count --miRNA ${miRNA_name} --biotype 'promoter_plus' --output ${miRNA_name}_promoter_plus_annotation.txt

    """
}
