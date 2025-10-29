#!/usr/bin/env python3

import argparse
import tqdm
import sys
import pysam


def process(psvs, vcf, padding):
    SENTINEL = object()

    regions = []
    for psv in tqdm.tqdm(psvs):
        positions = [(psv.chrom, psv.start - padding, psv.start + len(psv.ref) + padding)]
        pos2 = psv.info.get('pos2')
        if pos2 is None:
            continue
        for entry in pos2:
            entry = entry.split(':')
            assert len(entry) <= 4
            allele_ix = int(entry[3]) if len(entry) == 4 else 1
            pos = int(entry[1]) - 1
            positions.append((entry[0], pos - padding, pos + len(psv.alleles[allele_ix]) + padding))
        assert len(positions) > 1

        trivial = True
        for chrom, start, end in positions:
            try:
                it = vcf.fetch(chrom, start, end)
            except ValueError:
                # Invalid contig
                continue
            if next(it, SENTINEL) is not SENTINEL:
                trivial = False
                break

        if not trivial:
            regions.extend(positions)

    contigs = { contig: i for i, contig in enumerate(psvs.header.contigs) }
    def sorter(tup):
        chrom, start, end = tup
        if chrom in contigs:
            return (contigs[chrom], '', start, end)
        return (sys.maxsize, chrom, start, end)

    regions.sort(key=sorter)
    return regions


def write_out(regions, out):
    curr_chrom = None
    curr_start = None
    curr_end = None
    for chrom, start, end in regions:
        start = max(start, 0)
        if curr_chrom == chrom and curr_end >= start:
            curr_end = max(curr_end, end)
        else:
            if curr_chrom is not None:
                out.write(f'{curr_chrom}\t{curr_start}\t{curr_end}\n')
            curr_chrom = chrom
            curr_start = start
            curr_end = end
    if curr_chrom is not None:
        out.write(f'{curr_chrom}\t{curr_start}\t{curr_end}\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--psvs', metavar='FILE', required=True,
        help='VCF file with PSVs.')
    parser.add_argument('-v', '--vars', metavar='FILE', required=True,
        help='VCF file for comparison. Only presence of variants is checked.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output BED file with non-trivial PSV sites.')
    parser.add_argument('-P', '--padding', metavar='INT', type=int, default=0,
        help='Padding [%(default)s].')
    args = parser.parse_args()

    psvs = pysam.VariantFile(args.psvs)
    vcf = pysam.VariantFile(args.vars)
    regions = process(psvs, vcf, args.padding)
    out = open(args.output, 'w')
    write_out(regions, out)


if __name__ == '__main__':
    main()
