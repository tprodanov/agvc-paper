import math
from itertools import product

"""
CODE for NCF1 LD calculation, fixed ploidy for both SNPS
"""

def multinom_coeff(counts):
    n = sum(counts)
    c = math.factorial(n)
    for x in counts:
        c //= math.factorial(x)
    return c

def feasible_states(n, GA, GB):
    """Enumerate all (n00,n01,n10,n11) consistent with totals."""
    states = []
    for n00, n01, n10 in product(range(n+1), repeat=3):
        n11 = n - (n00 + n01 + n10)
        if n11 < 0:
            continue
        # allele count constraints
        if n10 + n11 != GA: continue
        if n01 + n11 != GB: continue
        states.append((n00, n01, n10, n11))
    return states

def ld_EM_unphased(GA_list, GB_list, n=6, tol=1e-8, max_iter=200):
    """
    GA_list, GB_list: lists of allele-1 counts for SNP A and B for each sample.
    n = ploidy (default 6).
    Returns hap frequencies + LD stats.
    """
    samples = len(GA_list)

    # Precompute feasible states for each sample
    all_states = []
    for GA, GB in zip(GA_list, GB_list):
        all_states.append(feasible_states(n, GA, GB))

    # Initialize haplotype frequencies
    f = [0.25, 0.25, 0.25, 0.25]  # f00,f01,f10,f11

    for it in range(max_iter):
        exp_counts = [0.0, 0.0, 0.0, 0.0]  # expected total over samples

        for states in all_states:
            # compute denominator
            weights = []
            for (n00, n01, n10, n11) in states:
                counts = (n00, n01, n10, n11)
                w = multinom_coeff(counts) * \
                    (f[0]**n00) * (f[1]**n01) * (f[2]**n10) * (f[3]**n11)
                weights.append(w)
            total_w = sum(weights)

            # accumulate expected counts
            for w, (n00, n01, n10, n11) in zip(weights, states):
                p = w / total_w
                exp_counts[0] += p * n00
                exp_counts[1] += p * n01
                exp_counts[2] += p * n10
                exp_counts[3] += p * n11

        # M-step: normalize
        total = sum(exp_counts)
        new_f = [c / total for c in exp_counts]

        # convergence
        if max(abs(new_f[i] - f[i]) for i in range(4)) < tol:
            f = new_f
            break
        f = new_f

    # Extract haplotype frequencies
    f00, f01, f10, f11 = f

    # Compute allele freqs
    pA = f10 + f11
    pB = f01 + f11

    # LD parameter D
    D = f11 - pA * pB

    # r^2
    denom = pA * (1 - pA) * pB * (1 - pB)
    r2 = (D * D) / denom if denom > 0 else None

    # D'
    if D == 0:
        Dprime = 0
    else:
        if D > 0:
            Dmax = min(pA*(1-pB), (1-pA)*pB)
        else:
            Dmax = min(pA*pB, (1-pA)*(1-pB))
        Dprime = D / Dmax if Dmax > 0 else None

    return {
        "f00": f00, "f01": f01, "f10": f10, "f11": f11,
        "pA": pA, "pB": pB,
        "D": D, "Dprime": Dprime, "r2": r2
    }

GA = []
GB = []
with open('R90H_deltaGT.genotypes') as F:
    for line in F: 
        gen = line.strip().split()
        GA.append(gen[0].count('1'))
        GB.append(gen[1].count('1'))
#print(GA,GB)
#GA = [0, 1, 6, 3]  # allele-1 counts at SNP A per sample
#GB = [0, 0, 6, 2]  # allele-1 counts at SNP B per sample

result = ld_EM_unphased(GA, GB, n=6)
print(result)

