#!/usr/bin/env python3

import os
import sys
import itertools
import collections
import operator
import argparse
import pysam
import gzip
import numpy as np
from intervaltree import IntervalTree

import parascopy.inner.common as common


def load_regions(filename):
    trees = collections.defaultdict(IntervalTree)
    regions = set()
    with common.open_possible_gzip(filename) as inp:
        for line in inp:
            line = line.strip().split('\t')
            region = line[3] if len(line) > 3 else '+'
            trees[line[0]].addi(int(line[1]), int(line[2]), region)
            regions.add(region)
    regions = sorted(regions)
    regions.append('*')
    return trees, regions


def get_partial_counts(chrom, start, end, all_trees):
    l = end - start
    overlaps = collections.defaultdict(lambda: np.zeros(l, dtype=bool))
    raw_overlaps = []
    for trees in all_trees:
        trees[]
    for overlap in tree.overlap(start, end):
        overlaps[overlap.data][max(overlap.begin - start, 0) : min(overlap.end - start, l)] = True

    uncovered = np.ones(l, dtype=bool)
    for locus, coverage in overlaps.items():
        uncovered &= ~coverage
        yield locus, np.sum(coverage) / l
    unc = np.sum(uncovered)
    if unc:
        yield '*', unc / l


def get_overlap_counts(start, end, tree):
    overlap = set(map(operator.attrgetter('data'), tree.overlap(start, end)))
    if overlap:
        return ((locus, 1) for locus in overlap)
    return (('*', 1),)


def count(vcf_filename, all_trees, count_homologous, partial):
    """
    For all covered tuple with regions,
        returns pair (total number, number of SNPs, number of indels) for any and high qualities.
    """
    covered_positions = collections.defaultdict(IntervalTree) if count_homologous else None
    counts = collections.defaultdict(lambda: [0] * 6)
    with pysam.VariantFile(vcf_filename) as vcf:
        for record in vcf:
            ref_len = len(record.ref)
            regions = [(record.chrom, record.start, record.start + ref_len)]
            if count_homologous:
                for reg in record.info['pos2']:
                    if reg != '???':
                        try:
                            reg = reg.split(':')
                            start2 = int(reg[1]) - 1
                            regions.append((reg[0], start2, start2 + ref_len))
                        except IndexError:
                            sys.stderr.write(
                                f'\n\n\nError in {vcf_filename} with {record.chrom}:{record.start + 1}: {reg}\n\n\n')
                            raise

            high_qual = record.samples[0].get('GQ', 10000) >= 10
            is_snp = all(len(alt) == ref_len for alt in record.alts)

            for chrom, start, end in regions:
                assert start < end
                if covered_positions is not None:
                    if covered_positions[chrom].overlaps(start, end):
                        # This variant was already covered in one of the other copies.
                        continue
                    covered_positions[chrom].addi(start, end, None)

                it = get_partial_counts(start, end, trees[chrom]) if partial and start + 1 < end \
                    else get_overlap_counts(start, end, trees[chrom])
                for locus, size in it:
                    counts[locus][0] += size
                    counts[locus][2 - int(is_snp)] += size
                    if high_qual:
                        counts[locus][3] += size
                        counts[locus][5 - int(is_snp)] += size
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--eval', metavar='DIR', required=True,
        help='RTG evaluation directory.')
    parser.add_argument('-R', '--regions', metavar='FILE', required=True, nargs='+',
        help='BED file(s), where fourth column declares region type.')
    parser.add_argument('-o', '--output', metavar='FILE', required=True,
        help='Output CSV file.')
    parser.add_argument('-a', '--all', action='store_true',
        help='Output counts for all entries, even empty.')
    parser.add_argument('--partial', action='store_true',
        help='Count variants partially if they overlap a region partially.')
    parser.add_argument('--homologous', action='store_true',
        help='Account for homologous coordinates (in pos2 info field).')
    args = parser.parse_args()

    all_trees = []
    all_regions = []
    for regions_file in args.regions:
        trees, regions = load_regions(args.regions)
        all_trees.append(trees)
        all_regions.append(regions)

    counts_tpb = count(os.path.join(args.eval, 'tp-baseline.vcf.gz'), all_trees, args.homologous, args.partial)
    counts_tpc = count(os.path.join(args.eval, 'tp.vcf.gz'), all_trees, args.homologous, args.partial)
    counts_fp = count(os.path.join(args.eval, 'fp.vcf.gz'), all_trees, args.homologous, args.partial)
    counts_fn = count(os.path.join(args.eval, 'fn.vcf.gz'), all_trees, args.homologous, args.partial)

    types = ('any\tall', 'any\tsnps', 'any\tindels', 'high\tall', 'high\tsnps', 'high\tindels')
    with common.open_possible_gzip(args.output, 'w') as out:
        out.write('# {}\n'.format(' '.join(sys.argv)))
        out.write('\t'.join(f'region{i+1}' for i in range(len(all_regions))))
        out.write('\tqual\tvar_type\ttp_base\ttp_call\tfp\tfn\n')

        for curr_regions in itertools.product(*all_regions):
            for i, ty in enumerate(types):
                tpb = counts_tpb[curr_regions][i]
                tpc = counts_tpc[curr_regions][i]
                fp = counts_fp[curr_regions][i]
                fn = counts_fn[curr_regions][i]
                if i == 0 and not args.all and tpb + tpc + fp + fn == 0:
                    break
                out.write(f'{"\t".join(curr_regions)}\t{ty}\t{tpb:.5g}\t{tpc:.5g}\t{fp:.5g}\t{fn:.5g}\n')


if __name__ == '__main__':
    main()
