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

def complementMap = [
    'A': 'U',
    'T': 'A',
    'G': 'C',
    'C': 'G'
]

workflow READ_MIRNA_TABLE{
    take:
    csv_path

    main:
    Channel
        .fromPath(csv_path)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                miRNA: row.miRNA,
                mature: row.mature,
                seed: row.mature[1..6],
                target: row.mature[1..6].collect { base -> complementMap[base] }.join().reverse(),
                extract_seq: row.mature[2..-1],
                id: "${row.miRNA}_${row.mature}",
            ]

            return meta
        }
        .set { mirna_ch }
    
    emit:
    mirna_ch
}
