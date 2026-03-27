#!/usr/bin/env python3
"""
gen_genome_frags.py — Simulate aDNA reads as random fragments of a reference genome.

Unlike gen_synthetic (which generates isolated molecules), this tool samples reads
as overlapping fragments from a shared reference, so longer reads provide k-mers
that extend beyond shorter reads' termini.

This is the correct model for testing fqdup extend:
  - Each read = fragment [start, start+length) of the reference
  - C→T damage applied at the 5' end (position p from terminus)
  - G→A damage applied at the 3' end
  - Extension works because reads from different positions overlap each other

Usage:
  python3 gen_genome_frags.py [options] > out.fq

Options:
  --genome-len N     Reference genome length (default: 2000)
  --n-reads N        Number of reads to generate (default: 50000)
  --min-len N        Min fragment length (default: 30)
  --max-len N        Max fragment length (default: 80)
  --dmax5 F          C→T amplitude at 5' end (default: 0.0)
  --dmax3 F          G→A amplitude at 3' end (default: 0.0)
  --lambda5 F        5' decay rate (default: 0.35)
  --lambda3 F        3' decay rate (default: 0.35)
  --seed N           RNG seed (default: 42)
  --qual N           Base quality score (default: 37)
"""

import sys
import math
import random
import argparse

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--genome-len', type=int, default=2000)
    p.add_argument('--n-reads',    type=int, default=50000)
    p.add_argument('--min-len',    type=int, default=30)
    p.add_argument('--max-len',    type=int, default=80)
    p.add_argument('--dmax5',      type=float, default=0.0)
    p.add_argument('--dmax3',      type=float, default=0.0)
    p.add_argument('--lambda5',    type=float, default=0.35)
    p.add_argument('--lambda3',    type=float, default=0.35)
    p.add_argument('--seed',       type=int, default=42)
    p.add_argument('--qual',       type=int, default=37)
    return p.parse_args()

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')

def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]

def apply_damage(seq, dmax5, dmax3, lam5, lam3, rng):
    """Apply C→T at 5' and G→A at 3' with exponential decay."""
    seq = list(seq)
    L = len(seq)
    for p in range(min(L, 15)):
        prob = dmax5 * math.exp(-lam5 * p)
        if seq[p] in ('C', 'c') and rng.random() < prob:
            seq[p] = 'T'
    for p in range(min(L, 15)):
        dist = L - 1 - p    # distance from 3' end
        if dist < 0:
            break
        prob = dmax3 * math.exp(-lam3 * p)
        if seq[dist] in ('G', 'g') and rng.random() < prob:
            seq[dist] = 'A'
    return ''.join(seq)

def main():
    args = parse_args()
    rng = random.Random(args.seed)

    # Generate reference genome
    bases = 'ACGT'
    genome = ''.join(rng.choice(bases) for _ in range(args.genome_len))

    sys.stderr.write(
        f"gen_genome_frags: genome_len={args.genome_len} n_reads={args.n_reads} "
        f"min_len={args.min_len} max_len={args.max_len} "
        f"dmax5={args.dmax5} dmax3={args.dmax3}\n"
    )

    qual_char = chr(args.qual + 33)
    out = sys.stdout

    for i in range(args.n_reads):
        # Sample fragment length and start position
        flen = rng.randint(args.min_len, args.max_len)
        start = rng.randint(0, args.genome_len - flen)
        seq = genome[start:start + flen]

        # Randomly orient (simulate both strands)
        if rng.random() < 0.5:
            seq = revcomp(seq)

        # Apply damage
        if args.dmax5 > 0 or args.dmax3 > 0:
            seq = apply_damage(seq, args.dmax5, args.dmax3,
                               args.lambda5, args.lambda3, rng)

        qual = qual_char * len(seq)
        out.write(f'@r{i+1}\n{seq}\n+\n{qual}\n')

if __name__ == '__main__':
    main()
