process BEDTOOLS_CALL_PEAK {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.filter.bed"), emit: bed

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def bin         = params.bin
    def new_meta    = meta + [single: true]
    def miRNA       = meta.miRNA
    def mismatch    = params.extract_mismatch
    def reads_cut   = params.cluster_reads_cut

    """
    # 1. bam to bed
    bedtools bamtobed -i ${bam} > ${prefix}.tmp.bed

    # 2. sort and cluster
    sort -k1,1 -k2,2n ${prefix}.tmp.bed | bedtools cluster -d 100 -i - > ${prefix}.clustered.bed

    # 3. merge and filter clusters 
    awk -v reads_cut=${reads_cut} '
    {
        cluster = \$NF
        strand = \$6
        key = \$1"\t"strand"\t"cluster
        if (!(key in min_start)) {
            min_start[key] = \$2
            max_end[key] = \$3
            count[key] = 1
        } else {
            if (\$2 < min_start[key]) min_start[key] = \$2
            if (\$3 > max_end[key]) max_end[key] = \$3
            count[key]++
        }
    }
    END {
        for (key in min_start) {
            if (count[key] >= reads_cut) {
                split(key, arr, "\t")
                chr = arr[1]
                strand = arr[2]
                cluster = arr[3]
                print chr, min_start[key], max_end[key], cluster, count[key], strand
            }
        }
    }' OFS="\t" ${prefix}.clustered.bed > ${prefix}.cluster_summary.${reads_cut}.filter.bed

    """
}
