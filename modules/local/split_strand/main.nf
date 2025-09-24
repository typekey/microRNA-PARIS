
process SPLIT_STRAND {
    tag "$meta.id"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/macs3:3.0.1--py311h0152c62_3':
    //     'biocontainers/macs3:3.0.1--py311h0152c62_3' }"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("*.minus.bam"), emit: minus_bam
    tuple val(meta), path("*.plus.bam"), emit: plus_bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandness = meta.strandness ? meta.strandness : 'R' // default to RF if not specified
    
    """
    if [[ "$strandness" == "R" || "$strandness" == "F" ]]; then
        # R/F library
        # forward read
        samtools view -hub -f 0 -F 16 "$bam_file" | samtools sort -T ./ -@ "${task.cpus}" -o "${prefix}.plus.bam"
        # reverse read
        samtools view -hub -f 16 "$bam_file" | samtools sort -T ./ -@ "${task.cpus}" -o "${prefix}.minus.bam"
    elif [[ "$strandness" == "RF" ]]; then
        # RF library
        echo "RF library"
        echo "second in pair, forward read"
        # second in pair, forward read
        samtools view -@ "${task.cpus}" -hub -f 128 -F 16 -o "${prefix}.plus.bam.r2" "$bam_file"
        echo "first in pair, reverse read"
        # first in pair, reverse read
        samtools view -@ "${task.cpus}" -hub -f 80 -o "${prefix}.plus.bam.r1" "$bam_file"
        echo "merge"
        samtools merge -f "${prefix}.plus.bam" "${prefix}.plus.bam.r1" "${prefix}.plus.bam.r2"
        rm "${prefix}.plus.bam.r1" "${prefix}.plus.bam.r2"

        # second in pair, reverse read
        echo "second in pair, reverse read"
        samtools view -@ "${task.cpus}" -hub -f 144 -o "${prefix}.minus.bam.r2" "$bam_file"
        # first in pair, forward read
        echo "first in pair, forward read"
        samtools view -@ "${task.cpus}" -hub -f 64 -F 16 -o "${prefix}.minus.bam.r1" "$bam_file"
        echo "merge"
        samtools merge -f "${prefix}.minus.bam" "${prefix}.minus.bam.r1" "${prefix}.minus.bam.r2"
        rm "${prefix}.minus.bam.r1" "${prefix}.minus.bam.r2"
    else
        # FR library
        # second in pair, reverse read
        samtools view -@ "${task.cpus}" -hub -f 144 -o "${prefix}.plus.bam.r1" "$bam_file"
        # first in pair, forward read
        samtools view -@ "${task.cpus}" -hub -f 64 -F 16 -o "${prefix}.plus.bam.r2" "$bam_file"
        samtools merge -f "${prefix}.plus.bam" "${prefix}.plus.bam.r1" "${prefix}.plus.bam.r2"
        rm "${prefix}.plus.bam.r1" "${prefix}.plus.bam.r2"

        # second in pair, forward read
        samtools view -@ "${task.cpus}" -hub -f 128 -F 16 -o "${prefix}.minus.bam.r1" "$bam_file"
        # first in pair, reverse read
        samtools view -@ "${task.cpus}" -hub -f 80 -o "${prefix}.minus.bam.r2" "$bam_file"
        samtools merge -f "${prefix}.minus.bam" "${prefix}.minus.bam.r1" "${prefix}.minus.bam.r2"
        rm "${prefix}.minus.bam.r1" "${prefix}.minus.bam.r2"
    fi
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
