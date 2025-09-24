
process MACS3_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macs3:3.0.1--py311h0152c62_3':
        'biocontainers/macs3:3.0.1--py311h0152c62_3' }"

    input:
    tuple val(meta), path(ipbam), val(dummy)
    val(strandedness)

    output:
    tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(meta), path("*.xls")                   , emit: xls
    path  "versions.yml"                             , emit: versions

    tuple val(meta), path("*.gappedPeak"), optional:true, emit: gapped
    tuple val(meta), path("*.bed")       , optional:true, emit: bed
    tuple val(meta), path("*.bdg")       , optional:true, emit: bdg

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()
    def format    = meta.single_end ? 'BAM' : 'BAMPE'
    def save_macs_pileup = params.save_macs_pileup ? '--bdg --SPMR' : ''
    def macs_pvalue = params.macs3_pvalue      ? "--pvalue ${params.macs3_pvalue}" : ''
    def macs_fdr = params.macs3_fdr         ? "--qvalue ${params.macs3_fdr}" : ''
    def macs3_gsize = params.effective_genome
    def keep_dup = params.keep_dup ? params.keep_dup.toString() : 'auto'
    def fragment_size = params.fragment_size ? "--nomodel --extsize ${params.fragment_size}" : ''
    def read_length = params.read_length ? "--tsize ${params.read_length}" : ''

    if(args_list.contains('--format')){
        def id = args_list.findIndexOf{it=='--format'}
        format = args_list[id+1]
        args_list.remove(id+1)
        args_list.remove(id)
    }
    """
    macs3 \\
        callpeak \\
        $save_macs_pileup \\
        $macs_pvalue \\
        $macs_fdr \\
        --keep-dup $keep_dup \\
        $fragment_size \\
        $read_length \\
        ${args_list.join(' ')} \\
        --gsize $macs3_gsize \\
        --name ${prefix}.${strandedness} \\
        --treatment $ipbam \\
        --nomodel \\
        --shift -25

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed -e "s/macs3 //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gappedPeak
    touch ${prefix}.bed
    touch ${prefix}.bdg
    touch ${prefix}.narrowPeak
    touch ${prefix}.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed -e "s/macs3 //g")
    END_VERSIONS
    """
}
