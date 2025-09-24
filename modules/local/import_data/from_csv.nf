#!/usr/bin/env nextflow
/*
----------------------------------------------------------------------------------------
File    :   from_csv.nf
Time    :   2025/01/07 23:11:54
Author  :   Lei Zheng
Version :   3.0
Contact :   type.zheng@gmail.com
Github  :   https://github.com/typekey
----------------------------------------------------------------------------------------
*/


workflow FROM_CSV{
    take:
    csv_path

    main:
    Channel
        .fromPath(csv_path)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                group: row.group,
                replicate: row.replicate,
                id: "${row.group}_${row.replicate}", // 创建唯一标识符
                single_end: false,  // single_end for miRNA chemirical reads
                strandedness: row.strandedness ? row.strandedness : 'R', // 默认值为 R
                miRNA: row.miRNA
            ]

            def reads = row.fastq_2 
                ? [file(row.fastq_1), file(row.fastq_2)] // PE 模式
                : [file(row.fastq_1)]                  // SE 模式
            
            return [meta, reads]
        }
        .set { fq_ch } // 输出 channel
    
    emit:
    data_ch = fq_ch
}

workflow FROM_CSV_PE{
    take:
    csv_path

    main:
    Channel
        .fromPath(csv_path)
        .splitCsv(header: true)
        .map {
            row ->
            tuple(row.group, row.replicate, file(row.fastq_1), file(row.fastq_2))
        }
        .set { fq_ch }
    
    emit:
    data_ch = fq_ch
}