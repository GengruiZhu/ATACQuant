#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATACQuant Phase 1 (Script 01): Topological Anchoring & Peak Clustering
Assigns empirical ATAC-seq peaks to GraphAllele syntenic allele clusters
via TSS proximity anchoring.

Ancestry labels in the input file are used as-is and must match the
--labels values supplied in downstream steps (no built-in translation).
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
        description="ATACQuant Topological Anchoring Engine",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-p", "--peaks",    required=True, help="MACS2 narrowPeak output file")
    parser.add_argument("-g", "--gff",      required=True, help="Genome annotation GFF3")
    parser.add_argument("-c", "--clusters", required=True, help="GraphAllele my_clusters.tsv")
    parser.add_argument("-a", "--ancestry", required=True, help="Gene ancestry mapping TSV (4th column = ancestry label)")
    parser.add_argument("--up",   type=int, default=3000, help="TSS upstream window (bp)")
    parser.add_argument("--down", type=int, default=1000, help="TSS downstream window (bp)")
    parser.add_argument("-o", "--outdir",   default="Topological_Clusters", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("=================================================================")
    print(f"[{get_now()}] ATACQuant: Topological Anchoring & Peak Clustering")

    # ================= 1. Load ancestry labels and allele clusters =================
    gene_ancestry = {}
    with open(args.ancestry, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('Gene_ID') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                # Use ancestry label as-is (uppercase); must match --labels downstream
                gene_ancestry[parts[0]] = parts[3].strip().upper()

    gene_to_cluster = {}
    with open(args.clusters, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            cid = parts[0]
            genes = parts[1].split(',')
            for g in genes:
                gene_to_cluster[g.strip()] = cid

    print(f"  -> Loaded {len(gene_ancestry)} gene ancestry labels, "
          f"{len(set(gene_to_cluster.values()))} allele clusters.")

    # ================= 2. Parse GFF3 and build TSS regulatory windows =================
    chr_genes = defaultdict(list)
    gene_count = 0

    with open(args.gff, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue

            chrom  = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            strand = parts[6]

            gene_id = None
            for item in parts[8].split(';'):
                if item.startswith('ID='):
                    gene_id = item.replace('ID=', '').strip()
                    break
            if not gene_id:
                continue

            if strand == '+':
                tss = start
                win_start = max(1, tss - args.up)
                win_end   = tss + args.down
            else:
                tss = end
                win_start = max(1, tss - args.down)
                win_end   = tss + args.up

            chr_genes[chrom].append((win_start, win_end, tss, gene_id))
            gene_count += 1

    print(f"  -> Extracted {gene_count} gene regulatory windows "
          f"(TSS -{args.up} ~ +{args.down} bp).")

    # Sort by window start for short-circuit scanning
    for chrom in chr_genes:
        chr_genes[chrom].sort(key=lambda x: x[0])

    # ================= 3. Anchor peaks to nearest gene within window =================
    em_clusters      = defaultdict(list)
    orphan_peaks     = []
    peak_ancestry_map = {}
    peak_coords_map  = {}

    peak_count  = 0
    mapped_count = 0

    with open(args.peaks, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts   = line.strip().split('\t')
            chrom   = parts[0]
            p_start = int(parts[1])
            p_end   = int(parts[2])
            peak_id = parts[3]
            peak_count += 1

            peak_coords_map[peak_id] = (chrom, p_start, p_end)

            p_mid     = (p_start + p_end) / 2.0
            best_gene = None
            min_dist  = float('inf')

            if chrom in chr_genes:
                for win_start, win_end, tss, gene_id in chr_genes[chrom]:
                    # Short-circuit: genes are sorted by win_start
                    if win_start > p_end:
                        break
                    if max(p_start, win_start) <= min(p_end, win_end):
                        dist = abs(p_mid - tss)
                        if dist < min_dist:
                            min_dist  = dist
                            best_gene = gene_id

            # ================= 4. Route peak to cluster or orphan list =================
            if best_gene:
                mapped_count += 1
                if best_gene in gene_to_cluster:
                    cid = gene_to_cluster[best_gene]
                    em_clusters[cid].append(peak_id)
                    peak_ancestry_map[peak_id] = gene_ancestry.get(best_gene, "UNKNOWN")
                else:
                    orphan_peaks.append((peak_id, chrom, p_start, p_end))
            else:
                orphan_peaks.append((peak_id, chrom, p_start, p_end))

    print(f"[{get_now()}] Routing complete. Total peaks: {peak_count} | Anchored: {mapped_count}")

    # ================= 5. Write output files =================
    cluster_out = os.path.join(args.outdir, "Syntenic_Peak_Clusters.tsv")
    saf_out     = os.path.join(args.outdir, "EM_Target_Peaks.saf")
    anc_out     = os.path.join(args.outdir, "Peak_Ancestry.txt")
    ls_out      = os.path.join(args.outdir, "Lineage_Specific_Peaks.tsv")

    valid_em_clusters = 0

    with open(cluster_out, 'w') as fc, open(saf_out, 'w') as fsaf, \
         open(anc_out, 'w') as fanc, open(ls_out, 'w') as fls:

        fsaf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        fanc.write("PeakID\tAncestry\n")

        # Write orphan peaks directly to lineage-specific file
        for op in orphan_peaks:
            fls.write(f"{op[0]}\t{op[1]}\t{op[2]}\t{op[3]}\n")

        # Clusters with >1 peak enter the EM deconvolution pool
        for cid, p_list in em_clusters.items():
            if len(p_list) > 1:
                fc.write(f"{cid}\t{','.join(p_list)}\n")
                valid_em_clusters += 1

                for pid in p_list:
                    chrom, start, end = peak_coords_map[pid]
                    fsaf.write(f"{pid}\t{chrom}\t{start}\t{end}\t.\n")
                    anc = peak_ancestry_map.get(pid, "UNKNOWN")
                    fanc.write(f"{pid}\t{anc}\n")
            else:
                # Single-peak cluster: no competition, route to lineage-specific
                pid = p_list[0]
                chrom, start, end = peak_coords_map[pid]
                fls.write(f"{pid}\t{chrom}\t{start}\t{end}\n")

    print("=================================================================")
    print(f"[SUCCESS] Topological anchoring complete.")
    print(f"  -> {valid_em_clusters} syntenic clusters queued for EM deconvolution.")
    print(f"  -> Lineage-specific peaks written to {ls_out}")
    print("=================================================================")

if __name__ == "__main__":
    main()
