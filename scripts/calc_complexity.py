#!/usr/bin/env python3

import argparse
import pysam
import sys


def calculate_complexity(seq, window, step, kmers):
    round_len = (len(seq) // step) * step
    divisors = [min(4 ** k, window - k + 1) for k in kmers]
    shift = (len(seq) - round_len) // 2

    for start in range(shift, len(seq) - window + 1, step):
        subseq = seq[start:start+window]
        assert len(subseq) == window
        complexity = 1.0
        for k, d in zip(kmers, divisors):
            complexity *= len(set(subseq[i:i+k] for i in range(0, window - k + 1))) / d
        yield start, complexity


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', metavar='FILE', required=True,
        help='Input reference fasta.')
    parser.add_argument('-R', '--regions', metavar='FILE', required=True,
        help='Either single region in format chrom:start-end or a BED file with regions.')
    parser.add_argument('-w', '--window', metavar='INT', default=100, type=int,
        help='Window size [%(default)s].')
    parser.add_argument('-s', '--step', metavar='INT', default=25, type=int,
        help='Step size [%(default)s].')
    parser.add_argument('-k', '--kmers', metavar='INT', nargs='+', default=[5], type=int,
        help='k-mer sizes, used for complexity calculation.')
    args = parser.parse_args()

    reference = pysam.FastaFile(args.reference)
    out = sys.stdout

    if ':' in args.regions:
        chrom, start_end = args.regions.split(':')
        start, end = start_end.split('-')
        f = [f'{chrom}\t{int(start)-1}\t{end}']
    else:
        f = open(args.regions)
    for line in f:
        if line.startswith('#'):
            continue
        chrom, start, end = line.strip().split('\t')[:3]
        start = int(start)
        end = int(end)
        seq = reference.fetch(chrom, start, end)
        for i, compl in calculate_complexity(seq, args.window, args.step, args.kmers):
            out.write(f'{chrom}\t{start + i}\t{start + i + args.window}\t{compl:.5f}\n')


if __name__ == '__main__':
    main()
