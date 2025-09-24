
process MERGE_STRAND {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(plus_peak), path(minus_peak)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.strand.narrowPeak"), emit: narrow_peak

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandness = meta.strandness
    
    """
    cat $plus_peak  | cut -f1-5 | awk -v OFS="\t" -v FS="\t" '{print \$0, "+"}'  > ${prefix}.bed
    cat $minus_peak | cut -f1-5 | awk -v OFS="\t" -v FS="\t" '{print \$0, "-"}' >> ${prefix}.bed

    sort -k1,1 -k2,2n "${prefix}.bed" -o "${prefix}.bed"

    awk -v OFS="\t" '{\$6 = "+"; print}' "$plus_peak" > "${prefix}.strand.narrowPeak"
    awk -v OFS="\t" '{\$6 = "-"; print}' "$minus_peak" >> "${prefix}.strand.narrowPeak"

    sort -k1,1 -k2,2n "${prefix}.strand.narrowPeak" -o "${prefix}.strand.narrowPeak"
    """
}
