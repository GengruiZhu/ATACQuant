#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATACQuant Phase 4 (Script 04): N-Way Bayesian EM Convergence Engine

Fuses three evidence streams:
  - K-mer prior (theta, from script 02b)
  - Per-read lineage likelihoods (P-values, from script 03)
  - Physical peak length (length-bias correction)

Runs Expectation-Maximization to deconvolve subgenome-specific
chromatin accessibility across syntenic homeologous peaks.

Supports N >= 2 subgenome species; species labels are inferred
automatically from the header of the P-value input file.
"""

import sys
import os
import argparse
from collections import defaultdict

def get_now():
    import time
    return time.strftime("%Y-%m-%d %H:%M:%S")

def main():
    parser = argparse.ArgumentParser(
        description="ATACQuant N-Way Bayesian EM Convergence Engine",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-c", "--clusters", required=True,
                        help="Phase 1 output: Syntenic_Peak_Clusters.tsv")
    parser.add_argument("-a", "--ancestry", required=True,
                        help="Phase 1 output: Peak_Ancestry.txt")
    parser.add_argument("-s", "--saf",      required=True,
                        help="Phase 1 output: EM_Target_Peaks.saf (for physical lengths)")
    parser.add_argument("-t", "--theta",    required=True,
                        help="Phase 2b output: Global_Theta_Matrix.tsv (k-mer prior)")
    parser.add_argument("-p", "--pvals",    required=True,
                        help="Phase 3 output: Reads_P_values.tsv (per-read likelihoods)")
    parser.add_argument("-o", "--out",    default="EM_Final_Expression.tsv",
                        help="Output file")
    parser.add_argument("--epsilon",  type=float, default=0.01,
                        help="Laplace smoothing pseudo-count")
    parser.add_argument("--max_iter", type=int,   default=300,
                        help="Maximum EM iterations")
    parser.add_argument("--tol",      type=float, default=1e-5,
                        help="Convergence tolerance")
    args = parser.parse_args()

    print("\n" + "=" * 65)
    print(f"[{get_now()}] ATACQuant: Bayesian EM Convergence Engine")

    # ================= 1. Load physical lengths and ancestry labels =================
    peak_lengths = {}
    with open(args.saf, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith(('GeneID', 'PeakID')) or not line.strip():
                continue
            parts = line.strip().split('\t')
            pid, start, end = parts[0], int(parts[2]), int(parts[3])
            peak_lengths[pid] = max(end - start + 1, 1)

    ancestry_dict = {}
    with open(args.ancestry, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith(('GeneID', 'PeakID')) or not line.strip():
                continue
            parts = line.strip().split('\t')
            ancestry_dict[parts[0]] = parts[1].upper()

    print(f"  -> Loaded {len(peak_lengths)} peak physical lengths and ancestry labels.")

    # ================= 2. Load cluster structure and theta priors =================
    cluster_alleles = defaultdict(list)
    peak_to_cluster = {}
    with open(args.clusters, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith(('#', 'Cluster')) or not line.strip():
                continue
            parts = line.strip().split('\t')
            cid   = parts[0]
            peaks = parts[1].split(',')
            for p in peaks:
                p_clean = p.strip()
                cluster_alleles[cid].append(p_clean)
                peak_to_cluster[p_clean] = cid

    initial_theta = {}
    with open(args.theta, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('Cluster') or not line.strip():
                continue
            parts = line.strip().split('\t')
            pid   = parts[1]
            if pid in peak_to_cluster:
                initial_theta[pid] = float(parts[5])

    print(f"  -> Anchored {len(cluster_alleles)} syntenic clusters with prior theta.")

    # ================= 3. Parse P-value file header to detect species labels =================
    pval_labels  = []
    with open(args.pvals, 'r', encoding='utf-8') as f:
        header_cols = f.readline().strip().split('\t')
        # Header format: Read_Name, Target_Peak, P_SP1, P_SP2, ..., Lineage_Status
        for col in header_cols[2:-1]:
            if col.startswith('P_'):
                pval_labels.append(col[2:])  # strip 'P_' prefix

    n_species     = len(pval_labels)
    label_to_idx  = {l: i for i, l in enumerate(pval_labels)}
    default_p     = 1.0 / n_species

    print(f"  -> Detected {n_species} subgenome species from P-value header: {pval_labels}")

    # ================= 4. Load P-values (memory-compressed into unique states) =================
    print(f"[{get_now()}] Loading per-read P-values into compressed state representation...")
    cluster_read_states = defaultdict(lambda: defaultdict(int))
    total_reads_loaded  = 0

    with open(args.pvals, 'r', encoding='utf-8') as f:
        f.readline()  # skip header
        for line in f:
            if not line.strip():
                continue
            parts       = line.strip().split('\t')
            target_peak = parts[1]

            if target_peak in peak_to_cluster:
                cid     = peak_to_cluster[target_peak]
                p_state = tuple(float(parts[2 + i]) for i in range(n_species))
                cluster_read_states[cid][p_state] += 1
                total_reads_loaded += 1

    print(f"  -> Compressed {total_reads_loaded:,} reads into unique probability states.")

    # ================= 5. EM iterations =================
    print(f"\n[{get_now()}] Launching Bayesian EM (length-corrected, {n_species}-way)...")
    final_results   = []
    total_clusters  = len(cluster_alleles)
    processed       = 0

    for cid, alleles in cluster_alleles.items():
        processed += 1
        if processed % 5000 == 0:
            print(f"  -> Deconvolved {processed}/{total_clusters} clusters "
                  f"({processed / total_clusters * 100:.1f}%)...")

        # --- Warm start: Laplace-smoothed theta from k-mer prior ---
        smoothed  = {a: (initial_theta.get(a, 0.0) + args.epsilon) for a in alleles}
        total_s   = sum(smoothed.values())
        theta     = {a: smoothed[a] / total_s for a in alleles}

        reads_info = cluster_read_states.get(cid, {})
        if not reads_info:
            for a in alleles:
                final_results.append((cid, a, ancestry_dict.get(a, "UNKNOWN"),
                                      peak_lengths.get(a, 1), 0.0, theta[a], theta[a]))
            continue

        # --- EM main loop ---
        for iteration in range(args.max_iter):
            old_theta       = theta.copy()
            expected_counts = {a: 0.0 for a in alleles}

            # E-Step: distribute reads using theta * length * lineage_prior
            for p_state, state_count in reads_info.items():
                weights    = {}
                weight_sum = 0.0

                for a in alleles:
                    lin   = ancestry_dict.get(a, "UNKNOWN")
                    idx   = label_to_idx.get(lin, -1)
                    prior_p = p_state[idx] if idx >= 0 else default_p

                    L_a      = peak_lengths.get(a, 1.0)
                    w        = theta[a] * L_a * prior_p
                    weights[a] = w
                    weight_sum += w

                if weight_sum > 0:
                    for a in alleles:
                        expected_counts[a] += (weights[a] / weight_sum) * state_count

            # M-Step: remove length bias, compute opening density
            density_sum = 0.0
            densities   = {}
            for a in alleles:
                L_a          = peak_lengths.get(a, 1.0)
                d            = (expected_counts[a] + args.epsilon) / L_a
                densities[a] = d
                density_sum += d

            for a in alleles:
                theta[a] = densities[a] / density_sum

            # Convergence check
            diff = sum(abs(theta[a] - old_theta[a]) for a in alleles)
            if diff < args.tol:
                break

        for a in alleles:
            final_results.append((cid, a, ancestry_dict.get(a, "UNKNOWN"),
                                  peak_lengths.get(a, 1),
                                  expected_counts[a], densities[a], theta[a]))

    # ================= 6. Write output =================
    with open(args.out, 'w', encoding='utf-8') as fout:
        fout.write("ClusterID\tPeakID\tAncestry\tLength_bp\t"
                   "EM_Allocated_Reads\tFinal_Density\tFinal_Theta\n")
        for res in final_results:
            fout.write(f"{res[0]}\t{res[1]}\t{res[2]}\t{res[3]}\t"
                       f"{res[4]:.2f}\t{res[5]:.6e}\t{res[6]:.6f}\n")

    print(f"  -> Deconvolved {processed}/{total_clusters} clusters (100.0%).")
    print(f"\n[SUCCESS] EM convergence complete.")
    print(f"  -> Length-corrected accessibility matrix written to: {args.out}")
    print("=" * 65 + "\n")

if __name__ == "__main__":
    main()
