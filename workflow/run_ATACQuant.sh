#!/bin/bash
# ==============================================================================
# ATACQuant cluster submission template
# Edit the parameters below and submit with: qsub run_ATACQuant.sh
# (or adapt the header for SLURM: #SBATCH directives)
# ==============================================================================
#PBS -N ATACQuant
#PBS -l nodes=1:ppn=48
#PBS -q comput
#PBS -o atacquant.log
#PBS -j oe

cd $PBS_O_WORKDIR

date -R
set -e

# ---------------------------------------------------------------------------- #
# INPUT FILES — edit these paths
# ---------------------------------------------------------------------------- #

R1=/path/to/clean_R1.fq
R2=/path/to/clean_R2.fq

FASTA=/path/to/genome.fasta
MASKED=/path/to/genome.fasta.masked
GFF=/path/to/genome.gff3

CLUSTERS=/path/to/my_clusters.tsv
MATRIX=/path/to/PolyAlleler_Global_Matrix_Cleaned.tsv
ANCESTRY=/path/to/Gene_Ancestry_Mapping.tsv

# One k-mer FASTA per subgenome species, matching --labels order
KMER_SP1=/path/to/Species1_k15.fasta
KMER_SP2=/path/to/Species2_k15.fasta
KMER_SP3=/path/to/Species3_k15.fasta   # remove or add lines as needed

# ---------------------------------------------------------------------------- #
# SETTINGS — edit as needed
# ---------------------------------------------------------------------------- #

SAMPLE=MySample    # used as output column name and file prefix
THREADS=48
KMER_K=15          # k-mer length; must match your k-mer FASTA files
OUTDIR=ATACQuant_Out

# ---------------------------------------------------------------------------- #
# RUN
# ---------------------------------------------------------------------------- #

python ATACQuant.py \
    --r1     "$R1" \
    --r2     "$R2" \
    -f       "$FASTA" \
    -m       "$MASKED" \
    -g       "$GFF" \
    -c       "$CLUSTERS" \
    --matrix "$MATRIX" \
    -a       "$ANCESTRY" \
    --kmers  "$KMER_SP1" "$KMER_SP2" "$KMER_SP3" \
    --labels SP1 SP2 SP3 \
    -t       "$THREADS" \
    -k       "$KMER_K" \
    --prefix "$SAMPLE" \
    -o       "$OUTDIR"

date -R
echo "ATACQuant finished."
