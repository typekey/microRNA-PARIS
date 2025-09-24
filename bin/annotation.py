#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
----------------------------------------------------------------------------------------
@File    :   annotation.py
@Time    :   2024/11/25 15:25:34
@Author  :   Lei Zheng
@Version :   1.0
@Contact :   type.zheng@gmail.com
----------------------------------------------------------------------------------------
'''

import pandas as pd
import argparse

def read_featurecount_file(file):
    count_df = pd.read_csv(file, sep='\t', comment='#')
    count_df.columns = [x.split('/')[-1].replace('.bam', '').replace('.sorted', '').replace('.markdup','') for x in count_df.columns]
    return count_df

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process feature count file and gene annotation file with miRNA sequence.')
    parser.add_argument('--count_file', type=str, required=True, help='Path to the feature count file')
    parser.add_argument('--gene_anno_file', type=str, default='/data/user/leizheng/reference/Homo_sapiens/GENCODE/GRCh38/annotation/gencode.v44.annotation.ercc.gene.tsv', help='Path to the gene annotation file')
    parser.add_argument('--miRNA', type=str, required=True, help='miRNA sequence')
    parser.add_argument('--output', type=str, default='./annotation.csv', help='Path of output file')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Print current parameters
    print(f"Feature count file: {args.count_file}")
    print(f"Gene annotation file: {args.gene_anno_file}")
    print(f"miRNA sequence: {args.miRNA}")
    
    # Load count and gene annotation data
    count_df = read_featurecount_file(args.count_file)
    gene_df = pd.read_csv(args.gene_anno_file, sep='\t')

    # Filter count data and merge with gene annotation
    count_filter_df = count_df[count_df[count_df.columns[-1]] > 0]
    res_df = pd.merge(count_filter_df, gene_df, left_on='Geneid', right_on='gene_id', how='inner').loc[:, ['Geneid', 'gene_name', 'gene_type', count_df.columns[-1]]]
    res_df['miRNA'] = args.miRNA
    res_df.columns = ['gene_id', 'gene_name', 'gene_type', 'counts', 'miRNA_seq']

    res_df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
