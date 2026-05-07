import pysam
import sys
import numpy as np
from itertools import combinations
import networkx as nx

def get_allele_counts(gt, ploidy):
    """Count number of alternate alleles."""
    if gt is None or None in gt:
        return None
    return sum(1 for allele in gt if allele == 1)

def calculate_maf(allele_counts, ploidies):
    """Calculate MAF for a variant."""
    total_alleles = np.nansum(ploidies)
    alt_alleles = np.nansum(allele_counts)
    if total_alleles == 0:
        return np.nan
    freq = alt_alleles / total_alleles
    return min(freq, 1 - freq)

def load_variant_genotypes(vcf_path, maf_threshold=0.01, chrom=None, region=None):
    """Load genotypes from a VCF and apply MAF filter."""
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    variants = []
    positions = []
    macs = []
    pos_index = {}
    v=0

    for rec in vcf.fetch(chrom, region[0], region[1]) if chrom and region else vcf.fetch():
        if not rec.alts or len(rec.alts) > 1:
            continue  # Skip multiallelic
        
        gts = []
        ploidies = []
        for sample in samples:
            call = rec.samples[sample]
            gt = call['GT']
            ploidy = len(gt) if gt is not None else 2
            ac = get_allele_counts(gt, ploidy)
            gts.append(ac if ac is not None else np.nan)
            ploidies.append(ploidy)
        
        gts = np.array(gts)
        ploidies = np.array(ploidies)
        maf = calculate_maf(gts, ploidies)
        if rec.pos == 70951946: 
            for s in range(len(samples)):
                if gts[s] == 1 and ploidies[s] ==3: print(samples[s])
            #print('var',rec.pos,maf,gts,ploidies)

        if np.isnan(maf) or maf < maf_threshold:
            continue

        mac = int(np.nansum(gts))
        variants.append(gts)
        positions.append(rec.pos)
        macs.append(round(maf,3))
        pos_index[rec.pos] = v
        v +=1
    
    return np.array(variants), positions, macs, pos_index


variants, positions, macs,pos_index = load_variant_genotypes(sys.argv[1], maf_threshold=0.05)
