process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.6--py39hf95cd2a_1' :
        'biocontainers/cutadapt:4.6--py39hf95cd2a_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim_miRNA.trim_adapterumi.fastq.gz'), emit: reads
    tuple val(meta), path('*.cutadapt.json')          , emit: json
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.trim.cutadapt.fastq.gz -p ${prefix}_2.trim.cutadapt.fastq.gz"
    def cutadapt_params = params.cutadapt_params ?: ''
    def miRNA = meta.miRNA

    """
    # Pre-alignment trimming and filtering
    cutadapt \
        --times 3 \
        -e 0.1 \
        -O 1 \
        --quality-cutoff 6 \
        -o ${prefix}.trim_miRNA.fastq \
        -g ${miRNA} \
        -j $task.cpus \
        --json ${prefix}.trim_miRNA.json \\
        $reads

    # Extract sequences from FASTQ
    # cat ${prefix}.trim_miRNA.fastq | awk 'NR%4==2' > ${prefix}.trim_miRNA.fastq.seq

    # Adapter trimming
    cutadapt \
        -O 1 \
        --times 3 \
        -e 0.1 \
        --quality-cutoff 6 \
        -m 18 \
        --discard-untrimmed \
        -o ${prefix}.trim_miRNA.trim_adapterumi.fastq \
        ${cutadapt_params} \
        --cut -1 \
        -j $task.cpus \
        --json ${prefix}.cutadapt.json \
        ${prefix}.trim_miRNA.fastq

    gzip ${prefix}.trim_miRNA.trim_adapterumi.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "${prefix}.trim.fastq.gz" : "${prefix}_1.trim.fastq.gz ${prefix}_2.trim.fastq.gz"
    """
    touch ${prefix}.cutadapt.log
    touch ${trimmed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}