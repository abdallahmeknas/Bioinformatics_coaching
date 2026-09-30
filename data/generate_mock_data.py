#!/usr/bin/env python3
"""
Generate a tiny, deterministic mock sequencing dataset for teaching.

Creates:
  reference.fa         a random 2 contig genome (50 kb total)
  truth_snps.tsv       the SNPs we planted, with each sample's true genotype
  <sample>_R1/R2.fastq.gz   simulated paired end Illumina style reads
  samplesheet.csv      nf-core style sheet: sample,fastq_1,fastq_2

Samples:
  sampleA  good quality
  sampleB  good quality
  sampleC  poor quality at read ends and more adapter read-through
           (great for showing what FastQC and fastp catch)

Usage:
  python data/generate_mock_data.py --outdir data/mock
  python data/generate_mock_data.py --outdir data/mock --coverage 50 --seed 7
"""

import argparse
import csv
import gzip
import random
from pathlib import Path

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

CONTIGS = [("chr1", 30_000), ("chr2", 20_000)]

# name, quality at the 3' end of the read, fraction of short fragments with adapter
SAMPLES = [
    ("sampleA", 34, 0.02),
    ("sampleB", 33, 0.02),
    ("sampleC", 14, 0.15),
]


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def make_reference(rng):
    return {name: "".join(rng.choice(BASES) for _ in range(length)) for name, length in CONTIGS}


def plant_snps(rng, reference, n_snps, margin=300):
    """Pick SNP sites and give every sample a genotype (0/0, 0/1 or 1/1)."""
    sites = []
    per_contig = {name: n_snps * len(seq) // sum(len(s) for s in reference.values())
                  for name, seq in reference.items()}
    for name, seq in reference.items():
        positions = sorted(rng.sample(range(margin, len(seq) - margin), per_contig[name]))
        for pos in positions:
            ref_base = seq[pos]
            alt_base = rng.choice([b for b in BASES if b != ref_base])
            genotypes = {}
            for sample, _, _ in SAMPLES:
                genotypes[sample] = rng.choices(["0/0", "0/1", "1/1"], weights=[0.3, 0.4, 0.3])[0]
            # Make sure every planted site is variant in at least one sample
            if all(gt == "0/0" for gt in genotypes.values()):
                genotypes[rng.choice([s for s, _, _ in SAMPLES])] = "0/1"
            sites.append((name, pos, ref_base, alt_base, genotypes))
    return sites


def build_haplotypes(reference, sites, sample):
    """Return two haplotypes (dict contig -> sequence) for a diploid sample."""
    haps = [{n: list(s) for n, s in reference.items()} for _ in range(2)]
    for contig, pos, _ref, alt, genotypes in sites:
        gt = genotypes[sample]
        if gt == "0/1":
            haps[1][contig][pos] = alt
        elif gt == "1/1":
            haps[0][contig][pos] = alt
            haps[1][contig][pos] = alt
    return [{n: "".join(s) for n, s in h.items()} for h in haps]


def add_errors(rng, seq, q_start, q_end):
    """Quality declines linearly along the read; errors follow the qualities."""
    length = len(seq)
    out_seq, out_qual = [], []
    for i, base in enumerate(seq):
        q = q_start + (q_end - q_start) * i / max(1, length - 1)
        q = int(max(2, min(40, q + rng.gauss(0, 2))))
        if rng.random() < 10 ** (-q / 10):
            base = rng.choice([b for b in BASES if b != base])
        out_seq.append(base)
        out_qual.append(chr(q + 33))
    return "".join(out_seq), "".join(out_qual)


def simulate_sample(rng, haplotypes, sample, q_end, adapter_rate, n_pairs,
                    read_len, frag_mean, frag_sd, outdir):
    contig_names = list(haplotypes[0].keys())
    weights = [len(haplotypes[0][n]) for n in contig_names]
    r1_path = outdir / f"{sample}_R1.fastq.gz"
    r2_path = outdir / f"{sample}_R2.fastq.gz"

    with gzip.open(r1_path, "wt") as r1_out, gzip.open(r2_path, "wt") as r2_out:
        for i in range(1, n_pairs + 1):
            hap = haplotypes[rng.randrange(2)]
            contig = rng.choices(contig_names, weights=weights)[0]
            seq = hap[contig]

            if rng.random() < adapter_rate:
                frag_len = rng.randint(60, read_len - 5)  # short insert: reads run into adapter
            else:
                frag_len = int(max(read_len + 20, rng.gauss(frag_mean, frag_sd)))

            start = rng.randrange(0, len(seq) - frag_len)
            fragment = seq[start:start + frag_len]
            if rng.random() < 0.5:
                fragment = revcomp(fragment)

            read1 = (fragment + ADAPTER_R1)[:read_len]
            read2 = (revcomp(fragment) + ADAPTER_R2)[:read_len]
            # pad very short inserts with random bases after the adapter
            read1 += "".join(rng.choice(BASES) for _ in range(read_len - len(read1)))
            read2 += "".join(rng.choice(BASES) for _ in range(read_len - len(read2)))

            s1, q1 = add_errors(rng, read1, 37, q_end)
            s2, q2 = add_errors(rng, read2, 36, q_end - 2)
            name = f"{sample}_{i:06d}"
            r1_out.write(f"@{name} 1:N:0:1\n{s1}\n+\n{q1}\n")
            r2_out.write(f"@{name} 2:N:0:1\n{s2}\n+\n{q2}\n")

    return r1_path.resolve(), r2_path.resolve()


def write_fasta(reference, path, width=60):
    with open(path, "w") as fh:
        for name, seq in reference.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="data/mock", help="output folder (default: data/mock)")
    parser.add_argument("--seed", type=int, default=42, help="random seed (default: 42)")
    parser.add_argument("--coverage", type=int, default=30, help="mean depth per sample (default: 30)")
    parser.add_argument("--read-length", type=int, default=100, help="read length (default: 100)")
    parser.add_argument("--snps", type=int, default=40, help="number of planted SNPs (default: 40)")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    reference = make_reference(rng)
    write_fasta(reference, outdir / "reference.fa")

    sites = plant_snps(rng, reference, args.snps)
    with open(outdir / "truth_snps.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\t" + "\t".join(s for s, _, _ in SAMPLES) + "\n")
        for contig, pos, ref, alt, gts in sites:
            fh.write(f"{contig}\t{pos + 1}\t{ref}\t{alt}\t" + "\t".join(gts[s] for s, _, _ in SAMPLES) + "\n")

    genome_len = sum(len(s) for s in reference.values())
    n_pairs = args.coverage * genome_len // (2 * args.read_length)

    rows = []
    for sample, q_end, adapter_rate in SAMPLES:
        haps = build_haplotypes(reference, sites, sample)
        r1, r2 = simulate_sample(rng, haps, sample, q_end, adapter_rate, n_pairs,
                                 args.read_length, 300, 30, outdir)
        rows.append({"sample": sample, "fastq_1": str(r1), "fastq_2": str(r2)})
        print(f"  {sample}: {n_pairs:,} read pairs")

    with open(outdir / "samplesheet.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Mock data written to {outdir}/ ({genome_len:,} bp genome, {len(sites)} planted SNPs)")


if __name__ == "__main__":
    main()
