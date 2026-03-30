#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATACQuant Phase 2 (Script 02): K-mer Dual-Layer Immune Dictionary Builder

Builds a genome-unique, transposon-free k-mer dictionary for each syntenic peak:
  Layer 1 — Hard-masked FASTA (RepeatMasker N-bases eliminate transposon k-mers)
  Layer 2 — Bowtie2 reverse alignment against the unmasked genome
             (any k-mer mapping >= 2 times is discarded as non-unique)

Output: a binary (.pkl) dictionary mapping each surviving k-mer to its
        (cluster_id, peak_id) of origin, plus per-peak probe counts.
"""

import sys
import os
import argparse
import subprocess
import pickle
from collections import defaultdict

try:
    from Bio import SeqIO
except ImportError:
    print("[FATAL ERROR] BioPython not found. Please run: pip install biopython")
    sys.exit(1)

def get_now():
    import time
    return time.strftime("%Y-%m-%d %H:%M:%S")

def extract_peak_sequences(saf_file, masked_fasta, out_fa):
    """Extract peak sequences from hard-masked genome using bedtools."""
    print(f"[{get_now()}] Extracting peak sequences from masked genome via bedtools...")

    try:
        subprocess.run(["bedtools", "--version"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except Exception:
        print("[FATAL ERROR] bedtools not found. Please ensure it is installed.")
        sys.exit(1)

    temp_bed = out_fa + ".tmp.bed"
    with open(saf_file, 'r') as fin, open(temp_bed, 'w') as fout:
        next(fin)  # skip header
        for line in fin:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            fout.write(f"{parts[1]}\t{parts[2]}\t{parts[3]}\t{parts[0]}\n")

    cmd = ["bedtools", "getfasta", "-fi", masked_fasta,
           "-bed", temp_bed, "-nameOnly", "-fo", out_fa]
    subprocess.run(cmd, check=True)
    os.remove(temp_bed)
    print(f"  -> Masked peak sequences extracted to {out_fa}")

def main():
    parser = argparse.ArgumentParser(
        description="ATACQuant K-mer Dual-Layer Immune Dictionary Builder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-s", "--saf",          required=True,
                        help="Phase 1 output: EM_Target_Peaks.saf")
    parser.add_argument("-m", "--masked_fasta", required=True,
                        help="Hard-masked genome FASTA (RepeatMasker output)")
    parser.add_argument("-x", "--index",        required=True,
                        help="Bowtie2 index prefix for the UNMASKED genome")
    parser.add_argument("-c", "--clusters",     required=True,
                        help="Phase 1 output: Syntenic_Peak_Clusters.tsv")
    parser.add_argument("-k", "--kmer",   type=int, default=31,
                        help="K-mer length (recommend 31 for large genomes)")
    parser.add_argument("-t", "--threads", type=int, default=16,
                        help="Bowtie2 parallel threads")
    parser.add_argument("-o", "--outdir", default="Immune_Dictionary",
                        help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    k = args.kmer

    print("=================================================================")
    print(f"[{get_now()}] ATACQuant: K-mer Dual-Layer Immune Dictionary Builder")
    print(f"  K-mer length : {k}")
    print(f"  Threads      : {args.threads}")

    # ================= 1. Load cluster-to-peak mapping =================
    peak_to_cluster = {}
    with open(args.clusters, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            cid   = parts[0]
            peaks = parts[1].split(',')
            for p in peaks:
                peak_to_cluster[p.strip()] = cid

    # ================= 2. Layer 1 — transposon masking via N-bases =================
    raw_peak_fa = os.path.join(args.outdir, "Target_Peaks.fa")
    extract_peak_sequences(args.saf, args.masked_fasta, raw_peak_fa)

    print(f"[{get_now()}] Fragmenting sequences into {k}-mers (N-base filter active)...")
    kmer_to_peaks = defaultdict(set)

    for record in SeqIO.parse(raw_peak_fa, "fasta"):
        peak_id = str(record.id).split('::')[0].split(':')[0]
        if peak_id not in peak_to_cluster:
            continue

        seq_str = str(record.seq).upper()
        if len(seq_str) < k:
            continue

        for i in range(len(seq_str) - k + 1):
            kmer = seq_str[i:i+k]
            # Layer 1: any k-mer containing N (masked transposon) is discarded
            if 'N' not in kmer:
                kmer_to_peaks[kmer].add(peak_id)

    candidate_kmers = {}
    for kmer, p_set in kmer_to_peaks.items():
        if len(p_set) == 1:
            candidate_kmers[kmer] = list(p_set)[0]

    del kmer_to_peaks  # free memory

    print(f"  -> Layer 1 complete: {len(candidate_kmers)} candidate unique k-mers retained.")

    # ================= 3. Layer 2 — genome uniqueness via Bowtie2 =================
    print(f"[{get_now()}] Running Bowtie2 reverse scan for genome-wide uniqueness check...")

    candidate_fa = os.path.join(args.outdir, "Candidate_Kmers.fa")
    with open(candidate_fa, 'w') as f:
        for kmer in candidate_kmers.keys():
            # Use k-mer string as sequence ID to avoid reverse-complement ambiguity
            f.write(f">{kmer}\n{kmer}\n")

    sam_out = os.path.join(args.outdir, "Kmers_Mapped.sam")

    cmd_bowtie2 = [
        "bowtie2", "-p", str(args.threads), "-x", args.index,
        "-f", candidate_fa, "-k", "2", "-L", str(k), "-N", "0",
        "--no-unal", "-S", sam_out
    ]
    subprocess.run(cmd_bowtie2, check=True)

    # ================= 4. Parse SAM: keep only singly-mapping k-mers =================
    print(f"[{get_now()}] Parsing alignment results and applying uniqueness filter...")

    kmer_hit_counts = defaultdict(int)
    with open(sam_out, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            parts = line.split('\t')
            kmer_hit_counts[parts[0]] += 1

    master_dict  = {}
    probe_counts = defaultdict(int)

    for kmer, pid in candidate_kmers.items():
        if kmer_hit_counts.get(kmer, 0) == 1:
            cid = peak_to_cluster[pid]
            master_dict[kmer] = (cid, pid)
            probe_counts[pid] += 1

    os.remove(candidate_fa)
    os.remove(sam_out)

    # ================= 5. Serialize dictionary =================
    dict_pkl = os.path.join(args.outdir, "ATAC_Master_Immune_Dict.pkl")
    with open(dict_pkl, 'wb') as f:
        pickle.dump((master_dict, probe_counts), f)

    print("=================================================================")
    print(f"[SUCCESS] Immune dictionary construction complete.")
    print(f"  -> {len(master_dict)} genome-unique, transposon-free k-mer probes retained.")
    print(f"  -> Dictionary serialized to: {dict_pkl}")
    print("=================================================================")

if __name__ == "__main__":
    main()
