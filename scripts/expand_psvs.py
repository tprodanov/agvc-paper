#!/usr/bin/env python3

import sys
import pysam
import argparse
import numpy as np

import parascopy.inner.common as common
import parascopy.inner.genome as genome_


def process_psvs(vcf, out_header):
    for var in vcf:
        pos2 = var.info.get('pos2')
        if pos2 is None:
            continue
        if 'fval' in var.info:
            fvals = list(map(float, var.info['fval']))
        else:
            fvals = np.full(len(pos2) + 1, np.nan)
        main_pos = '{}:{}'.format(var.chrom, var.pos)

        rec = out_header.new_record()
        rec.chrom = var.chrom
        rec.start = var.start
        rec.alleles = var.alleles
        rec.info['fval'] = fvals[0]
        rec.info['main'] = main_pos
        yield rec

        for i, pos in enumerate(pos2):
            pos = pos.split(':')
            strand = pos[2] == '+'
            ref_ix = 1 if len(pos) == 3 else int(pos[3])

            rec = out_header.new_record()
            rec.chrom = pos[0]
            rec.start = int(pos[1]) - 1
            rec.ref = common.cond_rev_comp(var.alleles[ref_ix], strand=strand)
            rec.alts = [common.cond_rev_comp(seq, strand=strand)
                for j, seq in enumerate(var.alleles) if j != ref_ix]
            rec.info['fval'] = fvals[i + 1]
            rec.info['main'] = main_pos
            yield rec


def create_header(genome):
    vcf_header = pysam.VariantHeader()
    vcf_header.add_line('##command="{}"'.format(' '.join(sys.argv)))

    for name, length in zip(genome.chrom_names, genome.chrom_lengths):
        vcf_header.add_line('##contig=<ID={},length={}>'.format(name, length))
    vcf_header.add_line('##INFO=<ID=fval,Number=1,Type=Float,Description="Reference allele frequency">')
    vcf_header.add_line('##INFO=<ID=main,Number=1,Type=String,Description="Main position">')
    return vcf_header


# def filter_overlaps(records):
#     count = 0
#     curr_start_ix = 0
#     curr_chrom = None
#     curr_end = None

#     for i, rec in enumerate(records):
#         if curr_chrom != rec.chrom or (curr_chrom == rec.chrom and curr_end > rec.start):
#             if i - 1 > curr_start_ix:
#                 count += i - curr_start_ix
#                 for j in range(curr_start_ix, i):
#                     records[j].info['fval'] = np.nan
#             curr_chrom = rec.chrom
#             curr_start_ix = i
#             curr_end = rec.start + len(rec.ref)
#         else:
#             curr_end = max(curr_end, rec.start + len(rec.ref))

#     if i - 1 > curr_start_ix:
#         count += i - curr_start_ix
#         for j in range(curr_start_ix, i):
#             records[j].info['fval'] = np.nan
#     sys.stderr.write('Filtered {:,} overlapping variants\n'.format(count))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs='+', required=True,
        help='Multiple input VCF files with PSVs.')
    parser.add_argument('-f', '--fasta-ref', required=True,
        help='Fasta reference file.')
    parser.add_argument('-o', '--output', required=True,
        help='Combined output VCF file.')
    args = parser.parse_args()

    genome = genome_.Genome(args.fasta_ref)
    header = create_header(genome)

    records = []
    for filename in args.input:
        records.extend(process_psvs(pysam.VariantFile(filename), header))
    records.sort(key=lambda record: (genome.chrom_id(record.chrom), record.start))
    # filter_overlaps(records)

    with pysam.VariantFile(args.output, 'wz' if args.output.endswith('.gz') else 'w', header=header) as out_vcf:
        for rec in records:
            out_vcf.write(rec)

    if args.output.endswith('.gz'):
        common.Process(['tabix', '-p', 'vcf', args.output]).finish()


if __name__ == '__main__':
    main()
