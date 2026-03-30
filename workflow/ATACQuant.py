#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATACQuant End-to-End Master Controller

Orchestrates the full ATACQuant pipeline across 7 sequential phases:
  00 — Alignment, deduplication, Tn5 shift, peak calling
  01 — Topological peak anchoring to GraphAllele clusters
  02 — K-mer immune dictionary construction (dual-layer)
  02b — Prior theta estimation from FASTQ k-mer scanning
  03 — Per-read N-way lineage P-value assignment
  04 — Bayesian EM deconvolution
  05 — Final count matrix assembly
  06 — Local gene annotation via GraphAllele matrix

Supports an arbitrary number of subgenome species (N >= 2) via --kmers / --labels.
Checkpoint-based resumption: re-run with --force to restart from scratch.
"""

import os
import sys
import argparse
import subprocess
import time

# Dynamically derive paths relative to this script
WORKFLOW_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(WORKFLOW_DIR)
BIN_DIR      = os.path.join(PROJECT_ROOT, "bin")

def get_now():
    return time.strftime("%Y-%m-%d %H:%M:%S")

def run_cmd(cmd_str, step_name, checkpoint_file):
    print(f"\n{'=' * 75}")
    print(f"[{get_now()}] >>> Starting: {step_name}")
    print(f"Command: {cmd_str}")
    try:
        subprocess.run(cmd_str, shell=True, check=True, executable='/bin/bash')
        with open(checkpoint_file, 'w') as f:
            f.write(f"SUCCESS at {get_now()}\n")
        print(f"[{get_now()}] <<< Completed: {step_name}")
    except subprocess.CalledProcessError:
        print(f"\n[FATAL ERROR] {step_name} failed. Pipeline aborted.")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        prog="ATACQuant.py",
        description="ATACQuant: Subgenome-resolved ATAC-seq quantification for polyploids",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Input data
    grp_data = parser.add_argument_group("Core input data")
    grp_data.add_argument("--r1",         required=True, metavar="R1.fq",
                          help="Clean paired-end FASTQ Read 1 (plain or gzipped)")
    grp_data.add_argument("--r2",         required=True, metavar="R2.fq",
                          help="Clean paired-end FASTQ Read 2 (plain or gzipped)")
    grp_data.add_argument("-f", "--fasta", required=True, metavar="genome.fasta",
                          help="Reference genome FASTA (unmasked)")
    grp_data.add_argument("-m", "--masked", required=True, metavar="genome.masked.fasta",
                          help="Hard-masked reference genome FASTA (RepeatMasker output)")
    grp_data.add_argument("-g", "--gff",   required=True, metavar="genome.gff3",
                          help="Genome annotation GFF3")

    # Topology and ancestry
    grp_topo = parser.add_argument_group("Polyploid topology and lineage knowledge base")
    grp_topo.add_argument("-c", "--clusters", required=True,
                          help="GraphAllele my_clusters.tsv")
    grp_topo.add_argument("--matrix",         required=True,
                          help="GraphAllele PolyAlleler_Global_Matrix_Cleaned.tsv")
    grp_topo.add_argument("-a", "--ancestry",  required=True,
                          help="Gene ancestry mapping TSV (4th column = ancestry label, "
                               "uppercase, must match --labels)")
    grp_topo.add_argument("--kmers",  nargs="+", required=True,
                          metavar="SP.fasta",
                          help="Species-specific diagnostic k-mer FASTA files "
                               "(one per subgenome, order must match --labels)")
    grp_topo.add_argument("--labels", nargs="+", required=True,
                          metavar="LABEL",
                          help="Short species labels matching --kmers order "
                               "(e.g. SP1 SP2 SP3); must match ancestry file labels")

    # TSS window (forwarded to script 01)
    grp_tss = parser.add_argument_group("TSS regulatory window (script 01)")
    grp_tss.add_argument("--up",   type=int, default=3000,
                         help="Upstream window from TSS (bp)")
    grp_tss.add_argument("--down", type=int, default=1000,
                         help="Downstream window from TSS (bp)")

    # System settings
    grp_sys = parser.add_argument_group("System and output settings")
    grp_sys.add_argument("-k", "--kmer",   type=int, default=15,
                         help="K-mer length for lineage P-value scoring (script 03); "
                              "must match your k-mer FASTA files")
    grp_sys.add_argument("-t", "--threads", default="48",
                         help="CPU threads for Bowtie2 / samtools")
    grp_sys.add_argument("--prefix", default="sample",
                         help="Sample name used in output filenames and matrix column header")
    grp_sys.add_argument("-o", "--outdir",
                         default=os.path.join(WORKFLOW_DIR, "ATACQuant_Out"),
                         help="Output directory")
    grp_sys.add_argument("--genome-prefix", default=None,
                         help="Genome name prefix in allele IDs for script 06 "
                              "(e.g. 'GENOME' for 'GENOME_Chr10A0001'); auto-detected if omitted")
    grp_sys.add_argument("--force", action="store_true",
                         help="Ignore checkpoints and re-run all steps from scratch")

    args = parser.parse_args()

    # ================= Validate inputs =================
    if len(args.kmers) != len(args.labels):
        print(f"[FATAL ERROR] --kmers ({len(args.kmers)}) and "
              f"--labels ({len(args.labels)}) must have the same count.")
        sys.exit(1)
    if len(args.labels) < 2:
        print("[FATAL ERROR] At least 2 species are required (--labels).")
        sys.exit(1)

    # ================= Checkpoint setup =================
    os.makedirs(args.outdir, exist_ok=True)
    ckpt_dir = os.path.join(args.outdir, ".checkpoints")
    os.makedirs(ckpt_dir, exist_ok=True)

    if args.force:
        for fn in os.listdir(ckpt_dir):
            if fn.endswith(".done"):
                os.remove(os.path.join(ckpt_dir, fn))

    def ckpt(step_id):
        return os.path.join(ckpt_dir, f"{step_id}.done")

    def is_done(step_id):
        return os.path.exists(ckpt(step_id))

    # ================= Resolve script paths =================
    scripts = {
        "00":  os.path.join(BIN_DIR, "00_running_mapping.sh"),
        "01":  os.path.join(BIN_DIR, "01_peak_clustering.py"),
        "02":  os.path.join(BIN_DIR, "02_build_dict.py"),
        "02b": os.path.join(BIN_DIR, "02b_scan_theta.py"),
        "03":  os.path.join(BIN_DIR, "03_read_Pvalue.py"),
        "04":  os.path.join(BIN_DIR, "04_EM_final.py"),
        "05":  os.path.join(BIN_DIR, "05_Matrix_Integration.py"),
        "06":  os.path.join(BIN_DIR, "06_Peak_Annotation.py"),
    }

    # ================= Define intermediate paths =================
    narrow_peak = os.path.join(args.outdir, f"{args.prefix}_peaks.narrowPeak")
    dedup_bam   = os.path.join(args.outdir, f"{args.prefix}.dedup.bam")
    syn_clusters = os.path.join(args.outdir, "Syntenic_Peak_Clusters.tsv")
    saf_file    = os.path.join(args.outdir, "EM_Target_Peaks.saf")
    anc_file    = os.path.join(args.outdir, "Peak_Ancestry.txt")
    lineage_tsv = os.path.join(args.outdir, "Lineage_Specific_Peaks.tsv")
    dict_pkl    = os.path.join(args.outdir, "ATAC_Master_Immune_Dict.pkl")
    theta_tsv   = os.path.join(args.outdir, "Global_Theta_Matrix.tsv")
    pvals_tsv   = os.path.join(args.outdir, "Reads_P_values.tsv")
    em_final    = os.path.join(args.outdir, "EM_Final_Expression.tsv")
    final_count = os.path.join(args.outdir, f"{args.prefix}_Final_Count_Matrix.tsv")
    annot_out   = os.path.join(args.outdir, f"{args.prefix}_Local_Annotated_Matrix.tsv")

    bt2_idx = os.path.splitext(args.fasta)[0]

    # Build --kmers and --labels strings for subprocess calls
    kmers_str  = " ".join(args.kmers)
    labels_str = " ".join(args.labels)

    # ================= Execute pipeline =================

    # [00] Alignment + peak calling
    if not is_done("00_mapping"):
        cmd = (f"bash {scripts['00']} "
               f"--r1 {args.r1} --r2 {args.r2} "
               f"--fasta {args.fasta} --prefix {args.prefix} "
               f"--threads {args.threads} --outdir {args.outdir}")
        run_cmd(cmd, "00 — Alignment & Peak Calling", ckpt("00_mapping"))

    # [01] Topological peak anchoring
    if not is_done("01_topo"):
        cmd = (f"python {scripts['01']} "
               f"-p {narrow_peak} -g {args.gff} "
               f"-c {args.clusters} -a {args.ancestry} "
               f"--up {args.up} --down {args.down} "
               f"-o {args.outdir}")
        run_cmd(cmd, "01 — Topological Peak Anchoring", ckpt("01_topo"))

    # [02] K-mer immune dictionary (locked at k=31 for genome-wide uniqueness)
    if not is_done("02_dict"):
        cmd = (f"python {scripts['02']} "
               f"-s {saf_file} -m {args.masked} -x {bt2_idx} "
               f"-c {syn_clusters} -k 31 -t {args.threads} -o {args.outdir}")
        run_cmd(cmd, "02 — K-mer Immune Dictionary", ckpt("02_dict"))

    # [02b] Prior theta from FASTQ (locked at k=31 to match dictionary)
    if not is_done("02b_theta"):
        cmd = (f"python {scripts['02b']} "
               f"-d {dict_pkl} -c {syn_clusters} "
               f"-r1 {args.r1} -r2 {args.r2} "
               f"-k 31 -o {theta_tsv}")
        run_cmd(cmd, "02b — Prior Theta Estimation", ckpt("02b_theta"))

    # [03] Per-read N-way lineage P-values (uses user-specified k for species k-mers)
    if not is_done("03_pvalue"):
        cmd = (f"python {scripts['03']} "
               f"-s {saf_file} -b {dedup_bam} "
               f"--kmers {kmers_str} --labels {labels_str} "
               f"-k {args.kmer} -o {pvals_tsv}")
        run_cmd(cmd, "03 — Lineage P-value Assignment", ckpt("03_pvalue"))

    # [04] Bayesian EM deconvolution
    if not is_done("04_em"):
        cmd = (f"python {scripts['04']} "
               f"-c {syn_clusters} -a {anc_file} -s {saf_file} "
               f"-t {theta_tsv} -p {pvals_tsv} -o {em_final}")
        run_cmd(cmd, "04 — Bayesian EM Deconvolution", ckpt("04_em"))

    # [05] Count matrix assembly
    if not is_done("05_matrix"):
        cmd = (f"python {scripts['05']} "
               f"-e {em_final} -l {lineage_tsv} -s {saf_file} "
               f"-b {dedup_bam} -n {args.prefix} -o {final_count}")
        run_cmd(cmd, "05 — Count Matrix Assembly", ckpt("05_matrix"))

    # [06] Local gene annotation
    if not is_done("06_annot"):
        gp_flag = f"--genome-prefix {args.genome_prefix}" if args.genome_prefix else ""
        cmd = (f"python {scripts['06']} "
               f"-i {final_count} -m {args.matrix} -y {syn_clusters} "
               f"-o {annot_out} {gp_flag}")
        run_cmd(cmd, "06 — Local Gene Annotation", ckpt("06_annot"))

    print(f"\n{'=' * 75}")
    print(f"[MISSION COMPLETE] ATACQuant pipeline finished successfully.")
    print(f"  Final annotated matrix : {annot_out}")
    print(f"  Timestamp              : {get_now()}")
    print("=" * 75)

if __name__ == "__main__":
    main()
