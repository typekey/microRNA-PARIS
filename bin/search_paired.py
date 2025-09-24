#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
----------------------------------------------------------------------------------------
@File    :   search_paired.py
@Time    :   2025/08/17 12:16:02
@Author  :   Lei Zheng
@Version :   1.0
@Contact :   type.zheng@gmail.com
----------------------------------------------------------------------------------------
'''

import os
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

def get_alignment(input_fastq, output_alignment, seed_sequence, max_mismatches=1, prefix_length=7, filtered_fastq=None):
    base_mapping = {"A": "U", "U": "A", "C": "G", "G": "C"}
    seed_sequence = seed_sequence.upper().replace("T", "U")
    prefix_seq = "".join([base_mapping[base] for base in seed_sequence[:prefix_length][::-1]])
    seed_length = len(seed_sequence)

    filtered_records = []
    seen_ids = set()

    with open(input_fastq) as handle, open(output_alignment, "w") as output:
        content = ''
        for record in SeqIO.parse(handle, "fastq"):
            input_sequence = str(record.seq).replace("T", "U")
            matched = False
            for i in range(len(input_sequence) - prefix_length + 1):
                window_seq = input_sequence[i:i + prefix_length]
                dist = distance(window_seq, prefix_seq)
                if dist <= max_mismatches:
                    matched = True
                    if i <= seed_length:
                        input_na = ' '*(seed_length-prefix_length-i)
                        middle_na = ' '*(seed_length-prefix_length)
                        seed_na = ' '*(i-seed_length+prefix_length)
                    else:
                        input_na = ''
                        middle_na = ' '*(i)
                        seed_na = ' '*(i-seed_length+prefix_length)
                    content += f"@{record.id}\n"
                    content += f"Target RNA : {input_na}{input_sequence}\n"
                    content += f"             {middle_na}{build_middle(input_sequence[i:i+prefix_length],seed_sequence[:prefix_length][::-1])}\n"
                    content += f"miRNA      : {seed_na}{seed_sequence[::-1]}\n\n"
                    break
            if matched and filtered_fastq and record.id not in seen_ids:
                filtered_records.append(record)
                seen_ids.add(record.id)
        output.write(content)
    if filtered_fastq and filtered_records:
        with open(filtered_fastq, "w") as fh:
            SeqIO.write(filtered_records, fh, "fastq")

def main():
    parser = argparse.ArgumentParser(description="miRNA seed alignment script")
    parser.add_argument("-i", "--input_fastq", required=True, help="Input FASTQ file")
    parser.add_argument("-o", "--output_alignment", required=True, help="Output alignment file")
    parser.add_argument("-s", "--seed_sequence", required=True, help="miRNA seed sequence")
    parser.add_argument("-m", "--max_mismatches", type=int, default=1, help="Maximum mismatches allowed")
    parser.add_argument("-p", "--prefix_length", type=int, default=7, help="Seed prefix length")
    parser.add_argument("-f", "--filtered_fastq", help="Output filtered FASTQ file")
    args = parser.parse_args()

    get_alignment(args.input_fastq, args.output_alignment, args.seed_sequence, args.max_mismatches, args.prefix_length, args.filtered_fastq)

if __name__ == "__main__":
    main()