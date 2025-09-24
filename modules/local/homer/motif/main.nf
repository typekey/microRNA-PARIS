
process HOMER_MOTIF {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bed_file)
    val(gene_bed)
    val(chrom_sizes)

    output:
    tuple val(meta), path("*.log")                   , emit: log
    tuple val(meta), path("*_homer")                 , emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gtf = params.gtf
    def fasta = params.fasta
    def homer_motif_length = params.homer_motif_length

    """
    sort -k5,5 -g ${bed_file} | awk 'FNR <= 2000{ print \$1"\\t"\$2"\\t"\$3}' > ${prefix}.location
    intersectBed -wo -a ${prefix}.location -b $gtf | awk -v OFS="\\t" '{print \$1,\$2,\$3,"*","*",\$10}' | sort -k1,2 | uniq > ${prefix}_bestpeaks.bed
    fastaFromBed -name+ -split -s -fi $fasta -bed ${prefix}_bestpeaks.bed > ${prefix}_bestpeaks.fa
    # ame -oc ${prefix}_ame ${prefix}_bestpeaks.fa m6A_motif.meme
    shuffleBed -incl ${gene_bed} -seed 12345 -noOverlapping -i ${prefix}_bestpeaks.bed -g ${chrom_sizes} > ${prefix}_random_peak.bed
    fastaFromBed -name+ -split -s -fi $fasta -bed ${prefix}_random_peak.bed > ${prefix}_random_peak.fa
    findMotifs.pl ${prefix}_bestpeaks.fa fasta ${prefix}_homer -fasta ${prefix}_random_peak.fa -p ${task.cpus} \
        -len $homer_motif_length -S 10 -rna -dumpFasta > ${prefix}_homer_run.log 2>&1
    """
}
