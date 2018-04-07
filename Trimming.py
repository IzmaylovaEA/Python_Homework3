
# coding: utf-8



#!/usr/bin/python
import argparse
import sys
from Bio import SeqIO
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='K-mers trimming')
    parser.add_argument('-i', '--input', help='Input FASTQ file', metavar='FILE', required=True)
    parser.add_argument('-t', '--threads', help='Number of threads', metavar='Int', type=int, default=1)
    parser.add_argument('-o', '--output', help='The resulting filename', metavar='FILE')
    parser.add_argument('-s', '--start', help = 'Number of cropped starting nucleotides', metavar='Int', type=int)
    parser.add_argument('-e', '--end', help = 'Number of cropped nucleotides at the end of reads', metavar='Int', type=int)
    parser.add_argument('-w', '--window', help = 'Sliding window', metavar='Int', type=int)
    parser.add_argument('-q', '--quality', help = 'Required mean quality', metavar='Int', type=int)
    args = parser.parse_args()
    sys.stdout = open(args.output, 'w')
    for record in SeqIO.parse(args.input, 'fastq'):
        seq_lng = len(record.seq)
        for index in range(args.start, seq_lng-args.end+1):
            s = 0
            end = 0
            for qual in record.letter_annotations['phred_quality'][index:(index+args.window)]:
                s += qual
            mean_qual = s/args.window
            if mean_qual < args.quality:
                end = index
                break
            else:
                end = index+args.window
        trimmed_read = record[args.start:end]
        SeqIO.write(trimmed_read, sys.stdout, 'fastq')
    sys.stdout.close()
