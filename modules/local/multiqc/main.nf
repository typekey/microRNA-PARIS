process MULTIQC {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.26--pyhdfd78af_0' :
        'biocontainers/multiqc:1.26--pyhdfd78af_0' }"

    input:
    path multiqc_config
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path logo
    path cutadapt_json
    path umi_extract_log
    path umi_dedup_log
    path ('fastqc/raw/*')
    path ('fastqc/trim/*')
    path ('trim_log/*')
    path ('star/*')
    path ('samtools/stats/*')
    path ('samtools/flagstat/*')
    path ('samtools/idxstats/*')
    path ch_fastq_screen_multiqc
    path ch_markduplicates_multiqc
    // path ch_frip_multiqc
    // path peak_count_multiqc

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    prefix = task.ext.prefix ?: "multiqc_report"
    """
    multiqc \\
        -n ${prefix}.html \\
        -f \\
        $args \\
        $custom_config \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
