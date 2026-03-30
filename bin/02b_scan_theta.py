#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATACQuant Phase 2b (Script 02b): ATAC-seq FASTQ Theta Scanner

Loads the immune k-mer dictionary produced by script 02 and scans both
FASTQ reads at high speed to compute the prior opening probability (Theta)
for each syntenic peak — used as the EM warm-start prior in script 04.
"""

import os
import sys
import gzip
import pickle
import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(
        description="ATACQuant K-mer FASTQ Theta Scanner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-d", "--dict",     required=True,
                        help="Immune dictionary produced by script 02 (ATAC_Master_Immune_Dict.pkl)")
    parser.add_argument("-c", "--clusters", required=True,
                        help="Phase 1 output: Syntenic_Peak_Clusters.tsv")
    parser.add_argument("-r1", "--read1",   required=True, help="Clean R1 FASTQ (plain or gzipped)")
    parser.add_argument("-r2", "--read2",   required=True, help="Clean R2 FASTQ (plain or gzipped)")
    parser.add_argument("-k", "--kmer",  type=int, default=31, help="K-mer length (must match script 02)")
    parser.add_argument("-o", "--out", default="Global_Theta_Matrix.tsv", help="Output file")
    args = parser.parse_args()

    # 1. Load immune dictionary
    with open(args.dict, 'rb') as f:
        master_dict, probe_counts = pickle.load(f)
    print(f"[INFO] Loaded {len(master_dict)} genome-unique immune probes.")

    # 2. Load cluster structure
    cluster_dict = defaultdict(list)
    with open(args.clusters, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            cluster_dict[parts[0]] = [p.strip() for p in parts[1].split(',')]

    # 3. Scan both FASTQ files (forward + reverse complement)
    hit_counts = defaultdict(int)
    trans = str.maketrans('ATCGN', 'TAGCN')

    def scan_fastq(fq_file):
        print(f"[INFO] Scanning {os.path.basename(fq_file)}...")
        open_fn = gzip.open if fq_file.endswith('.gz') else open

        # Pull globals into local scope to reduce per-iteration overhead
        k          = args.kmer
        local_dict = master_dict
        local_hits = hit_counts

        with open_fn(fq_file, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    seq = line.strip().upper()
                    seq_len = len(seq)
                    if seq_len < k:
                        continue

                    rc_seq = seq.translate(trans)[::-1]

                    for j in range(seq_len - k + 1):
                        kmer = seq[j:j+k]
                        if kmer in local_dict:
                            local_hits[local_dict[kmer][1]] += 1

                        rc_kmer = rc_seq[j:j+k]
                        if rc_kmer in local_dict:
                            local_hits[local_dict[rc_kmer][1]] += 1

                if i > 0 and i % 4_000_000 == 0:
                    print(f"  -> Processed {i // 4:,} reads...")

    scan_fastq(args.read1)
    scan_fastq(args.read2)

    # 4. Compute theta (prior opening probability) with zero-probability fix
    with open(args.out, "w") as out_f:
        out_f.write("Cluster_ID\tPeak_ID\tUnique_Probes\tRaw_Hits\tNormalized_Score\tTheta\n")

        for cid, peaks in cluster_dict.items():
            raw_scores  = {}
            total_score = 0.0

            for pid in peaks:
                probes = probe_counts.get(pid, 0)
                hits   = hit_counts.get(pid, 0)
                score  = (hits / probes) if probes > 0 else 0.0
                raw_scores[pid] = score
                total_score    += score

            for pid in peaks:
                # If all scores are zero, assign equal prior (avoids EM zero-probability deadlock)
                theta = (raw_scores[pid] / total_score) if total_score > 0 \
                        else (1.0 / len(peaks))
                out_f.write(
                    f"{cid}\t{pid}\t{probe_counts.get(pid, 0)}\t"
                    f"{hit_counts.get(pid, 0)}\t{raw_scores[pid]:.4f}\t{theta:.6f}\n"
                )

    print(f"[SUCCESS] Prior theta matrix written to: {args.out}")

if __name__ == "__main__":
    main()
