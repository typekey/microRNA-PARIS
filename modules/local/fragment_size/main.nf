
process MACS3_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(ipbam), path(controlbam)

    output:
    tuple val(meta), path("*fragment_size"), emit: fragment_size

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if(args_list.contains('--format')){
        def id = args_list.findIndexOf{it=='--format'}
        format = args_list[id+1]
        args_list.remove(id+1)
        args_list.remove(id)
    }
    """
    samtools view -@ ${task.cpus}  -f 2  ${ipbam} | awk '{if (\$9 > 0) print \$9}' | sort -n | awk '{a[i++]=\$1} END {print a[int(i/2)]}' > ${prefix}.fragment_size
    """
}
