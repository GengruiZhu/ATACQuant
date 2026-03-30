#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATACQuant Phase 3 (Script 03): N-Way Read Lineage Likelihood (P-value) Calculator

Scans the deduplicated BAM file at syntenic peak loci and assigns each
sequencing fragment a continuous Bayesian likelihood score for each subgenome
species, using species-specific diagnostic k-mer libraries.

Supports an arbitrary number of subgenomes (N >= 2).
Supply one k-mer FASTA per species via --kmers, and matching labels via --labels.
"""

import sys
import os
import argparse
import pysam

def get_now():
    import time
    return time.strftime("%Y-%m-%d %H:%M:%S")

def load_kmer_set(fasta_file, k):
    """Load a species-specific diagnostic k-mer set from a FASTA file."""
    kmer_set = set()
    if not os.path.exists(fasta_file):
        print(f"[FATAL ERROR] K-mer FASTA not found: {fasta_file}")
        sys.exit(1)

    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq = line.strip().upper()
                if len(seq) >= k:
                    kmer_set.add(seq[:k])
    return kmer_set

def calculate_read_P_nway(read_seq, kmer_sets, labels, k):
    """
    Assign lineage likelihood probabilities to a single read using N-species k-mer sets.
    Both forward and reverse-complement strands are scanned.

    Parameters
    ----------
    read_seq  : str
    kmer_sets : dict[label -> set]
    labels    : list[str]   ordered list of species labels
    k         : int         k-mer length

    Returns
    -------
    p_vals : dict[label -> float]   probability for each species (sum to ~1)
    status : str                    diagnostic string
    """
    n    = len(labels)
    hits = {label: 0 for label in labels}

    trans   = str.maketrans('ATCGN', 'TAGCN')
    fwd_seq = read_seq.upper()
    rev_seq = fwd_seq.translate(trans)[::-1]
    seq_len = len(fwd_seq)

    if seq_len < k:
        return {label: round(1.0 / n, 4) for label in labels}, "Too_Short"

    # Forward strand scan
    for i in range(seq_len - k + 1):
        kmer = fwd_seq[i:i+k]
        for label in labels:
            if kmer in kmer_sets[label]:
                hits[label] += 1

    # Reverse-complement scan
    for i in range(seq_len - k + 1):
        kmer = rev_seq[i:i+k]
        for label in labels:
            if kmer in kmer_sets[label]:
                hits[label] += 1

    total_hits = sum(hits.values())

    if total_hits > 0:
        eps         = 0.01
        denominator = total_hits + (n * eps)
        p_vals      = {label: (hits[label] + eps) / denominator for label in labels}

        max_label = max(hits, key=hits.get)
        max_val   = hits[max_label]
        n_max     = sum(1 for v in hits.values() if v == max_val)
        hit_str   = ":".join(str(hits[l]) for l in labels)

        if n_max == 1:
            status = f"{max_label}_Biased({hit_str})"
        else:
            status = f"Balanced/Tied({hit_str})"
    else:
        p_vals = {label: round(1.0 / n, 4) for label in labels}
        status = "Ambiguous_None"

    return p_vals, status

def main():
    parser = argparse.ArgumentParser(
        description="ATACQuant N-Way BAM Lineage P-value Assigner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-s", "--saf",    required=True,
                        help="Phase 1 output: EM_Target_Peaks.saf")
    parser.add_argument("-b", "--bam",    required=True,
                        help="Deduplicated BAM file (must retain query sequences)")
    parser.add_argument("--kmers",  nargs="+", required=True,
                        help="Species-specific k-mer FASTA files (one per species, order must match --labels)")
    parser.add_argument("--labels", nargs="+", required=True,
                        help="Species labels corresponding to --kmers (e.g. SP1 SP2 SP3)")
    parser.add_argument("-k", "--kmer",  type=int, default=21,
                        help="K-mer scan length (recommend 21; must match the k-mer FASTA files)")
    parser.add_argument("-o", "--out", default="Reads_P_values.tsv",
                        help="Output TSV file")
    args = parser.parse_args()

    # Validate label/kmer correspondence
    if len(args.kmers) != len(args.labels):
        print(f"[FATAL ERROR] --kmers ({len(args.kmers)} files) and "
              f"--labels ({len(args.labels)} labels) must have the same count.")
        sys.exit(1)
    if len(args.labels) < 2:
        print("[FATAL ERROR] At least 2 species are required (--labels).")
        sys.exit(1)

    labels = args.labels
    n      = len(labels)

    print("=================================================================")
    print(f"[{get_now()}] ATACQuant N-Way Lineage P-value Engine ({n} species: {', '.join(labels)})")

    # ================= 1. Load k-mer radar libraries =================
    print(f"[{get_now()}] Loading species-specific k-mer libraries...")
    kmer_sets = {}
    for label, fasta in zip(labels, args.kmers):
        kmer_sets[label] = load_kmer_set(fasta, args.kmer)
        print(f"  -> {label}: {len(kmer_sets[label]):,} probes loaded from {os.path.basename(fasta)}")

    # ================= 2. Parse SAF for target peak coordinates =================
    peak_coords = {}
    with open(args.saf, 'r') as f:
        for line in f:
            if line.startswith('GeneID') or not line.strip():
                continue
            parts   = line.strip().split('\t')
            peak_id = parts[0]
            peak_coords[peak_id] = (parts[1], int(parts[2]), int(parts[3]))

    print(f"  -> Locked {len(peak_coords)} syntenic target peaks.")

    # ================= 3. Scan BAM and assign P-values =================
    print(f"[{get_now()}] Scanning BAM file at syntenic peak loci...")

    try:
        bam_file = pysam.AlignmentFile(args.bam, "rb")
    except Exception as e:
        print(f"[FATAL ERROR] Cannot open BAM: {e}")
        sys.exit(1)

    total_reads_fetched  = 0
    total_reads_saved    = 0
    global_processed_reads = set()  # hash set for deduplication (memory-efficient)

    p_col_header = "\t".join(f"P_{l}" for l in labels)

    with open(args.out, "w") as out_f:
        out_f.write(f"Read_Name\tTarget_Peak\t{p_col_header}\tLineage_Status\n")

        for peak_id, (chrom, start, end) in peak_coords.items():
            try:
                fetch_iter = bam_file.fetch(chrom, max(0, start), end)
            except ValueError:
                continue

            for read in fetch_iter:
                # Only count Read1 to avoid double-counting fragments
                if not read.is_read1:
                    continue

                total_reads_fetched += 1

                if read.is_unmapped:
                    continue

                r_hash = hash(read.query_name)
                if r_hash in global_processed_reads:
                    continue

                seq = read.query_sequence
                if not seq:
                    continue

                global_processed_reads.add(r_hash)
                total_reads_saved += 1

                p_vals, status = calculate_read_P_nway(seq, kmer_sets, labels, args.kmer)
                p_str = "\t".join(f"{p_vals[l]:.4f}" for l in labels)

                out_f.write(f"{read.query_name}\t{peak_id}\t{p_str}\t{status}\n")

    bam_file.close()

    print("=================================================================")
    print(f"[REPORT] Total reads fetched (Read1 view) : {total_reads_fetched:,}")
    print(f"[REPORT] Unique fragments after global dedup: {total_reads_saved:,}")
    print(f"[SUCCESS] P-value matrix written to: {args.out}")
    print("=================================================================")

if __name__ == "__main__":
    main()
