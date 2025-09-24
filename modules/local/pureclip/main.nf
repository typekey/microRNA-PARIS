process PURECLIP {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}"

    input:
    tuple val(meta), path(bam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta = params.fasta
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools index $bam
    pureclip -i $bam -bai ${bam}.bai -g ${fasta} -iv 'chr1;chr2;chr3;' -nt 10 -o ${prefix}.crosslink_sites.bed
    """
}
