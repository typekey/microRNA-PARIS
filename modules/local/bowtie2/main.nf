#!/usr/bin/env nextflow
/*
----------------------------------------------------------------------------------------
File    :   bowtie2.nf
Time    :   2024/06/11 13:17:58
Author  :   Lei Zheng
Version :   1.0
Contact :   type.zheng@gmail.com
Github  :   https://github.com/typekey
----------------------------------------------------------------------------------------
*/


process BOWTIE2_MAPPING_GENOME {
    tag "${meta.id}"
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/genome/bowtie2", mode: 'copy'

    input:
    path fastq
    val meta
    
    output:
    val meta, emit: meta
    path '*.sorted.bam', emit: sorted_bam
    path '*.log', emit: log

    script:
    def bowtie2 = params.software.bowtie2
    def samtools = params.software.samtools
    def rmgenome_fq = meta.mode == "SE" ? "--un-gz ${meta.id}.rmgenome.fastq.gz" : "--un-conc-gz ${meta.id}.rmgenome.%.fastq.gz"
    def genome_bam = "${meta.id}.genome.align.sorted.bam"
    def index = params.reference[params.genome].bowtie2
    def input_fastq = meta.mode == "SE" ? "-U ${fastq[0]}" : "-1 ${fastq[0]} -2 ${fastq[1]}"
    def log = "${meta.id}.rmgenome.log"
    // --norc: the library type is fr-secondstrand.
    // --nofw (default): the library type is “fr-firststrand” and “unstranded”

    """
    ${bowtie2} -p ${task.cpus} -x ${index} \
    ${input_fastq} 2> ${log} \
    | ${samtools} view -@ ${task.cpus} -Shub /dev/stdin \
    | ${samtools} sort -T ${meta.id} -@ ${task.cpus} -o ${genome_bam} -
    """

    //  ${bowtie2} -p ${task.cpus} \
    // -X 2000 --local --no-mixed --no-discordant --no-unal \
    // -x ${index} \
    // ${input_fastq} ${rmgenome_fq} 2> ${log} \
    // | ${samtools} view -@ ${task.cpus} -Shub /dev/stdin \
    // | ${samtools} sort -T ${sample_id} -@ ${task.cpus} -o ${genome_bam} -
}

process BOWTIE2_MAPPING_CONTAMINATION {
    tag "${meta.id}"
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/contamination"

    input:
    path fastq
    val meta
    
    output:
    val meta, emit: meta
    path '*.fastq.gz', emit: rmcontamination_fastq
    // path '*.bam', emit: bam
    path '*.log', emit: log

    script:
    def bowtie2 = params.software.bowtie2
    def samtools = params.software.samtools
    def rmContamination_fq = meta.mode == "SE" ? "--un-gz ${meta.id}.rmcontamination.fastq.gz" : "--un-conc-gz ${meta.id}.rmcontamination.%.fastq.gz"
    def contamination_bam = "${meta.id}.contamination.align.sorted.bam"
    def index = params.reference[params.genome].contamination.bowtie2
    def input_fastq = meta.mode == "SE" ? "-U ${fastq[0]}" : "-1 ${fastq[0]} -2 ${fastq[1]}"
    def log = "${meta.id}.rmcontamination.log"
    def strict_paired_mode = params.remove_reads.strict_paired_mode ? "-X 2000 --local --no-mixed --no-discordant --no-unal":""
    // --norc: the library type is fr-secondstrand.
    // --nofw (default): the library type is “fr-firststrand” and “unstranded”

    """
    ${bowtie2} -p ${task.cpus} \
    ${strict_paired_mode} \
    -x ${index} \
    ${input_fastq} ${rmContamination_fq} -S ${meta.id}_contamination.sam 2> ${log}
    """

    // ${bowtie2} -p ${task.cpus} \
    // -X 2000 --local --no-mixed --no-discordant --no-unal \
    // -x ${index} \
    // ${input_fastq} ${rmContamination_fq} -S eg1.sam 2> ${log} \
    // | ${samtools} view -@ ${task.cpus} -Shub /dev/stdin \
    // | ${samtools} sort -T ${meta.id} -@ ${task.cpus} -o ${contamination_bam} -
}

process BOWTIE2_MAPPING_RRNA {
    tag "${meta.id}"
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/rRNA"

    input:
    path fastq
    val meta
    
    output:
    val meta, emit: meta
    path '*.fastq.gz', emit: rmrRNA_fastq
    path '*.sam', emit: sam
    path '*.log', emit: log

    script:
    def bowtie2 = params.software.bowtie2
    def samtools = params.software.samtools
    def rmrRNA_fq = meta.mode == "SE" ? "--un-gz ${meta.id}.rmrRNA.fastq.gz" : "--un-conc-gz ${meta.id}.rmrRNA.%.fastq.gz"

    def rRNA_bam = "${meta.id}.rRNA.align.sorted.bam"
    def index = params.reference[params.genome].rRNA.bowtie2
    def input_fastq = meta.mode == "SE" ? "-U ${fastq[0]}" : "-1 ${fastq[0]} -2 ${fastq[1]}"
    def log = "${meta.id}.rmrRNA.log"
    def strict_paired_mode = params.remove_reads.strict_paired_mode ? "-X 2000 --local --no-mixed --no-discordant --no-unal":""
    // --norc: the library type is fr-secondstrand.
    // --nofw (default): the library type is “fr-firststrand” and “unstranded”

    """
    ${bowtie2} -p ${task.cpus} \
    ${strict_paired_mode} \
    -x ${index} \
    ${input_fastq} ${rmrRNA_fq} -S ${meta.id}_rrna.sam 2> ${log}
    """

    // ${bowtie2} -p ${task.cpus} \
    // -X 2000 --local --no-mixed --no-discordant --no-unal \
    // -x ${index} \
    // ${input_fastq} ${rmrRNA_fq} 2> ${log} \
    // | ${samtools} view -@ ${task.cpus} -Shub /dev/stdin \
    // | ${samtools} sort -T ${meta.id} -@ ${task.cpus} -o ${rRNA_bam} -
}

process BOWTIE2_MAPPING_SRNA {
    tag "${meta.id}"
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/sRNA"

    input:
    path fastq
    val meta
    
    output:
    val meta, emit: meta
    path '*.fastq.gz', emit: rmsRNA_fastq
    path '*.bam', emit: bam
    path '*.log', emit: log

    script:
    def bowtie2 = params.software.bowtie2
    def samtools = params.software.samtools
    def rmsRNA_fq = meta.mode == "SE" ? "--un-gz ${meta.id}.rmsRNA.fastq.gz" : "--un-conc-gz ${meta.id}.rmsRNA.%.fastq.gz"

    def sRNA_bam = "${meta.id}.sRNA.align.sorted.bam"
    def index = params.reference[params.genome].sRNA.bowtie2
    def input_fastq = meta.mode == "SE" ? "-U ${fastq[0]}" : "-1 ${fastq[0]} -2 ${fastq[1]}"
    def log = "${meta.id}.rmsRNA.log"
    def strict_paired_mode = params.remove_reads.strict_paired_mode ? "-X 2000 --local --no-mixed --no-discordant --no-unal":""
    // --norc: the library type is fr-secondstrand.
    // --nofw (default): the library type is “fr-firststrand” and “unstranded”

    """
    ${bowtie2} -p ${task.cpus} \
    ${strict_paired_mode} \
    -x ${index} \
    ${input_fastq} ${rmsRNA_fq} 2> ${log} \
    | ${samtools} view -@ ${task.cpus} -Shub /dev/stdin \
    | ${samtools} sort -T ${meta.id} -@ ${task.cpus} -o ${sRNA_bam} -
    """
}

process BOWTIE2_MAPPING_SPIKE {
    tag "${meta.id}"
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/spike"

    input:
    path fastq_1
    path fastq_2

    output:
    path '*.1.fq', emit: fastq_1
    path '*.2.fq', emit: fastq_2
    path '*.sam', emit: sam
    path '*.report'

    script:
    { sampleId = fastq_1.toString().split("_")[0..1].join("_") }
    """
    ${params.tools.bowtie2} -p ${params.threads.bowtie2} \
    --nofw --no-unal --end-to-end -L 16 -N 1 --mp 5 \
    --un-conc ${sampleId}_spike.fq -x ${params.references.spike_expand.bowtie2} \
    -1 ${fastq_1} -2 ${fastq_2} > ${sampleId}_spike.sam 2> >(tee ${sampleId}_output.report >&2)
    """
}

process BOWTIE2_MAPPING_SNCRNA {
    tag "${meta.id}"
    debug params.debug
    label 'process_high'
    publishDir "${params.outdir.MAPPING}/sncRNA"

    input:
    path fastq_1
    path fastq_2

    output:
    path '*.1.fq', emit: fastq_1
    path '*.2.fq', emit: fastq_2
    path '*.sam', emit: sam
    path '*.report'

    script:
    { sampleId = fastq_1.toString().split("_")[0..1].join("_") }
    """
    ${params.tools.bowtie2} -p ${params.threads.bowtie2} \
    --nofw --all --no-unal --end-to-end -L 16 -N 1 --mp 5 \
    --un-conc ${sampleId}_sncRNA.fq -x ${params.references.sncRNA_human.bowtie2} \
    -1 ${fastq_1} -2 ${fastq_2} > ${sampleId}_sncRNA.sam 2> >(tee ${sampleId}_output.report >&2)
    """
}

process BOWTIE2_INDEX {
    tag "${meta.id}"
    debug params.debug
    label 'process_high'

    input:
    path fasta
    val index_name

    // output:
    // path '*'

    script:
    def index= index_name.toString()

    """
    {params.path_bowtie2} -p {params.threads} \
    --no-unal --end-to-end --fast \
    --un-conc {params.un} -x {params.ref_bowtie2} \
    -1 {input.r1} -2 {input.r2} > {output.sam} 2> >(tee {output.report} >&2)
    """
}