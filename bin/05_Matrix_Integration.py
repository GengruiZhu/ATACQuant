#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATACQuant Phase 5 (Script 05): Final Count Matrix Integration

Merges two peak categories into a single DESeq2-ready count matrix:
  1. Syntenic peaks — read counts from Bayesian EM deconvolution (script 04)
  2. Lineage-specific peaks — raw fragment counts from BAM (no deconvolution needed)

Physical genomic coordinates are attached to all peaks in the output.
"""

import os
import sys
import argparse
import pysam

def get_now():
    import time
    return time.strftime("%Y-%m-%d %H:%M:%S")

def main():
    parser = argparse.ArgumentParser(
        description="ATACQuant Final Count Matrix Builder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-e", "--em",      required=True,
                        help="Phase 4 output: EM_Final_Expression.tsv")
    parser.add_argument("-l", "--lineage", required=True,
                        help="Phase 1 output: Lineage_Specific_Peaks.tsv")
    parser.add_argument("-s", "--saf",     required=True,
                        help="Phase 1 output: EM_Target_Peaks.saf (for EM peak coordinates)")
    parser.add_argument("-b", "--bam",     required=True,
                        help="Deduplicated BAM file")
    parser.add_argument("-n", "--name",    default="Sample",
                        help="Sample column name in the output matrix")
    parser.add_argument("-o", "--out",     default="Final_Count_Matrix.tsv",
                        help="Output matrix file")
    args = parser.parse_args()

    print(f"[{get_now()}] ATACQuant: Final count matrix assembly...")

    # ================= 1. Load EM peak coordinates =================
    em_coords = {}
    with open(args.saf, 'r') as f:
        for line in f:
            if line.startswith('GeneID') or not line.strip():
                continue
            parts = line.strip().split('\t')
            em_coords[parts[0]] = (parts[1], parts[2], parts[3])  # Chr, Start, End

    print(f"  -> Loaded coordinates for {len(em_coords)} EM target peaks.")

    # ================= 2. Load EM deconvolved counts =================
    final_matrix = {}
    em_count = 0
    with open(args.em, 'r') as f:
        next(f)  # skip header
        for line in f:
            if not line.strip():
                continue
            parts    = line.strip().split('\t')
            peak_id  = parts[1]
            ancestry = parts[2]
            count    = int(round(float(parts[4])))
            chrom, start, end = em_coords.get(peak_id, ("Unknown", "0", "0"))
            final_matrix[peak_id] = (chrom, start, end, ancestry, "Syntenic_EM", count)
            em_count += 1

    print(f"  -> Loaded {em_count} EM-deconvolved syntenic peaks.")

    # ================= 3. Count raw fragments for lineage-specific peaks =================
    try:
        bam_file = pysam.AlignmentFile(args.bam, "rb")
    except Exception as e:
        print(f"[FATAL ERROR] Cannot open BAM: {e}")
        sys.exit(1)

    orphan_count = 0
    with open(args.lineage, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            parts   = line.strip().split('\t')
            peak_id = parts[0]
            chrom   = parts[1]
            start   = int(parts[2])
            end     = int(parts[3])
            ancestry = "UNKNOWN"

            try:
                fetch_iter     = bam_file.fetch(chrom, max(0, start), end)
                fragment_count = 0
                for read in fetch_iter:
                    # Mirror the EM engine's counting standard: Read1, mapped
                    if read.is_read1 and not read.is_unmapped:
                        fragment_count += 1

                final_matrix[peak_id] = (chrom, str(start), str(end),
                                         ancestry, "Orphan_Raw", fragment_count)
                orphan_count += 1
            except ValueError:
                continue

    bam_file.close()
    print(f"  -> Counted fragments for {orphan_count} lineage-specific peaks.")

    # ================= 4. Write merged output matrix =================
    with open(args.out, 'w') as fout:
        fout.write(f"PeakID\tChr\tStart\tEnd\tAncestry\tOrigin_Type\t{args.name}\n")
        for pid, (chrom, start, end, anc, origin, count) in final_matrix.items():
            fout.write(f"{pid}\t{chrom}\t{start}\t{end}\t{anc}\t{origin}\t{count}\n")

    print("=================================================================")
    print(f"[SUCCESS] Matrix assembly complete: {len(final_matrix)} regulatory elements.")
    print(f"  -> Output written to: {args.out}")
    print("=================================================================")

if __name__ == "__main__":
    main()
