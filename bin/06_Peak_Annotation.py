#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATACQuant Phase 6 (Script 06): Local Gene Precise Annotator

Maps each peak to its species-local allele gene ID using the GraphAllele
global matrix, by matching the peak's chromosome identity to the allele
entries within its syntenic cluster.

The optional --genome-prefix parameter handles genomes where chromosome IDs
are prefixed with a species name (e.g. 'GENOME_Chr10A' -> extracts 'Chr10A').
"""

import sys
import os
import argparse
import re

def get_now():
    import time
    return time.strftime("%Y-%m-%d %H:%M:%S")

def main():
    parser = argparse.ArgumentParser(
        description="ATACQuant Peak-to-Local-Gene Annotator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input",  required=True,
                        help="Phase 5 output: Final_Count_Matrix.tsv (must contain Chr column)")
    parser.add_argument("-m", "--matrix", required=True,
                        help="GraphAllele PolyAlleler_Global_Matrix_Cleaned.tsv")
    parser.add_argument("-y", "--syn",    required=True,
                        help="Phase 1 output: Syntenic_Peak_Clusters.tsv")
    parser.add_argument("-o", "--out",    required=True,
                        help="Output annotated matrix file")
    parser.add_argument("--genome-prefix", default=None,
                        help="Genome name prefix in allele IDs to handle "
                             "(e.g. 'GENOME' for IDs like 'GENOME_Chr10A0001'). "
                             "If omitted, any word prefix is stripped automatically.")
    args = parser.parse_args()

    print(f"[{get_now()}] ATACQuant Phase 06: Peak-to-Local-Gene Annotation")

    # Build chromosome extraction regex
    if args.genome_prefix:
        prefix_pat = re.escape(args.genome_prefix) + r'_'
    else:
        prefix_pat = r'\w+_'
    chr_pattern = re.compile(r'(?:' + prefix_pat + r')?(Chr\d+[A-Z]?)', re.IGNORECASE)

    # 1. Build PeakID -> ClusterID mapping
    peak_to_cluster = {}
    with open(args.syn, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            p   = line.strip().split('\t')
            cid = p[0]
            for pid in p[1].split(','):
                peak_to_cluster[pid.strip()] = cid

    print(f"  -> Loaded {len(peak_to_cluster)} peak-to-cluster mappings.")

    # 2. Build GraphAllele double index: cluster_id -> {core_chrom -> allele_id}
    cluster_index = {}
    with open(args.matrix, 'r') as f:
        f.readline()  # skip header
        for line in f:
            if not line.strip():
                continue
            p   = line.strip().split('\t')
            cid = p[0]
            cluster_index[cid] = {}
            for allele_id in p[3:]:  # columns Allele_A onwards
                if allele_id and allele_id != '-':
                    match = chr_pattern.search(allele_id)
                    if match:
                        core_chrom = match.group(1)
                        cluster_index[cid][core_chrom] = allele_id

    print(f"  -> Indexed {len(cluster_index)} GraphAllele syntenic clusters.")

    # 3. Annotate each peak in the count matrix
    mapped_genes = 0
    with open(args.input, 'r') as f, open(args.out, 'w') as fout:
        header_parts = f.readline().strip().split('\t')

        if header_parts[1] != "Chr":
            print("[FATAL ERROR] Column 2 of the input matrix is not 'Chr'. "
                  "Please use the output of Phase 5 without modification.")
            sys.exit(1)

        new_header = header_parts + ["ClusterID", "Local_Gene"]
        fout.write("\t".join(new_header) + "\n")

        for line in f:
            if not line.strip():
                continue
            parts      = line.strip().split('\t')
            pid        = parts[0]
            curr_chrom = parts[1]

            chrom_match    = chr_pattern.search(curr_chrom)
            core_curr_chrom = chrom_match.group(1) if chrom_match else curr_chrom

            curr_cluster = peak_to_cluster.get(pid, "Orphan")

            local_gene = "Intergenic/Distal"
            if curr_cluster in cluster_index:
                local_gene = cluster_index[curr_cluster].get(
                    core_curr_chrom, "No_Local_Gene_In_Cluster"
                )
                if local_gene not in ("Intergenic/Distal", "No_Local_Gene_In_Cluster"):
                    mapped_genes += 1

            fout.write(f"{line.strip()}\t{curr_cluster}\t{local_gene}\n")

    print(f"[{get_now()}] Annotation complete.")
    print(f"  -> {mapped_genes} peaks annotated with a local gene ID.")
    print(f"  -> Annotated matrix written to: {args.out}")
    print("=================================================================")

if __name__ == "__main__":
    main()
