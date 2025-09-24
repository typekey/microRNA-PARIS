#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
----------------------------------------------------------------------------------------
@File    :   build_mirna_paired.py
@Time    :   2025/03/18 11:20:12
@Author  :   Lei Zheng
@Version :   1.0
@Contact :   type.zheng@gmail.com
----------------------------------------------------------------------------------------
'''

import os
import gzip
import argparse
from Bio import SeqIO
from Levenshtein import distance

def build_middle(input_seq, seed):
    middle = ""
    base_mapping = {"A": "U", "U": "A", "C": "G", "G": "C"}
    for i in range(len(seed)):
        if input_seq[i] == base_mapping[seed[i]]:
            middle += "|"
        else:
            middle += "."
    return middle


def get_alignment(input_fastq, output_alignment, seed_sequence, max_mismatches = 1, prefix_length = 7):
    # base_mapping = {"A": "T", "T": "A", "C": "G", "G": "C"}
    base_mapping = {"A": "U", "U": "A", "C": "G", "G": "C"}
    seed_sequence = seed_sequence.upper().replace("T", "U")
    prefix_seq = "".join([base_mapping[base] for base in seed_sequence[:prefix_length][::-1]])
    seed_length = len(seed_sequence)
    output_motif_sequence = output_alignment.replace(".alignment", ".motif")

    with gzip.open(input_fastq, 'rt') as handle, gzip.open(output_alignment, "wt") as output, gzip.open(output_motif_sequence, "wt") as motif_output:
        content = ''
        motif_content=''
        content_simple = ''
        motif_seq=''

        # 判断是否为空
        for record in SeqIO.parse(handle, "fastq"):
            # print(str(record.seq))
            input_sequence = str(record.seq).replace("T", "U")
            input_length = len(input_sequence)
            for i in range(len(input_sequence) - prefix_length + 1):
                window_seq = input_sequence[i:i + prefix_length]
                dist = distance(window_seq, prefix_seq)
                if dist <= max_mismatches:
                    # print(i, prefix_seq, window_seq, dist)
                    if i <= seed_length:
                        input_na = ' '*(seed_length-prefix_length-i)
                        middle_na = ' '*(seed_length-prefix_length)
                        seed_na = ' '*(i-seed_length+prefix_length)
                        motif_seq = input_sequence[:i+prefix_length][-prefix_length:]
                        # motif_seq = motif_seq[:-prefix_length]
                        # print(">>>", motif_seq, motif_seq, prefix_length)
                    else:
                        input_na = ''
                        middle_na = ' '*(i)
                        seed_na = ' '*(i-seed_length+prefix_length)
                        motif_seq = input_sequence[i-seed_length:i+prefix_length][seed_length:]
                    content += f"@{record.id}\n"
                    content += f"Target RNA : {input_na}{input_sequence}\n"
                    content += f"             {middle_na}{build_middle(input_sequence[i:i+prefix_length],seed_sequence[:prefix_length][::-1])}\n"
                    content += f"miRNA      : {seed_na}{seed_sequence[::-1]}\n\n"
                    motif_content += motif_seq + "\n"
            
        output.write(content)
        motif_output.write(motif_content)


def main():
    parser = argparse.ArgumentParser(description="Align miRNA with target RNA sequences from a FASTQ file.")
    parser.add_argument("--input_fastq", help="Path to the input FASTQ file.")
    parser.add_argument("--output_dir", help="Directory to save the output alignment file.")
    parser.add_argument("--seed_sequence", help="Seed sequence for alignment.")
    parser.add_argument("--max_mismatches", type=int, default=1, help="Maximum number of mismatches allowed in the alignment.")
    parser.add_argument("--prefix_length", type=int, default=7, help="Length of the seed sequence prefix for alignment.")

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    output_alignment = os.path.join(args.output_dir, os.path.basename(args.input_fastq).replace(".fastq", ".alignment"))
    get_alignment(args.input_fastq, output_alignment, args.seed_sequence, args.max_mismatches, args.prefix_length)


if __name__ == "__main__":
    main()
