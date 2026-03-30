#!/bin/bash
# ==============================================================================
# ATACQuant Phase 00: Read Alignment, PCR Deduplication & Peak Calling
# Called by ATACQuant.py master controller; do not run directly.
# ==============================================================================
set -e

# ---------------------------------------------------------------------------- #
# Parse named arguments passed by ATACQuant.py                                 #
# ---------------------------------------------------------------------------- #
while [[ $# -gt 0 ]]; do
    case $1 in
        --r1)      R1="$2";      shift 2 ;;
        --r2)      R2="$2";      shift 2 ;;
        --fasta)   FASTA="$2";   shift 2 ;;
        --prefix)  PREFIX="$2";  shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --outdir)  OUT_DIR="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1"; exit 1 ;;
    esac
done

# Validate required arguments
for var in R1 R2 FASTA PREFIX THREADS OUT_DIR; do
    if [ -z "${!var}" ]; then
        echo "[FATAL ERROR] Missing required argument: --${var,,}"
        exit 1
    fi
done

mkdir -p "$OUT_DIR"

# Derive Bowtie2 index prefix from the FASTA path (strip extension)
BT2_IDX="${FASTA%.*}"

echo "=============================================================================="
echo "[Phase 00] ATACQuant: Alignment and Peak Calling"
echo "  R1           : $R1"
echo "  R2           : $R2"
echo "  Reference    : $FASTA"
echo "  Bowtie2 index: $BT2_IDX"
echo "  Prefix       : $PREFIX"
echo "  Threads      : $THREADS"
echo "  Output dir   : $OUT_DIR"
echo "=============================================================================="

# ================= 1. Build Bowtie2 index if absent =================
if [ ! -f "${BT2_IDX}.1.bt2" ] && [ ! -f "${BT2_IDX}.1.bt2l" ]; then
    echo "[STEP 1] Building Bowtie2 genome index..."
    bowtie2-build --threads "$THREADS" "$FASTA" "$BT2_IDX"
else
    echo "[STEP 1] Bowtie2 index already exists, skipping."
fi

# ================= 2. Alignment (retain multi-mappers for EM) =================
RAW_BAM="$OUT_DIR/${PREFIX}.raw.bam"
if [ ! -f "${RAW_BAM}.bai" ]; then
    echo "[STEP 2] Running Bowtie2 alignment (-k 20 retains full multi-mapping pool for EM)..."
    mkdir -p "$OUT_DIR/sort_tmp"

    bowtie2 -p "$THREADS" -x "$BT2_IDX" -1 "$R1" -2 "$R2" \
        -X 1000 --no-mixed --no-discordant -k 20 | \
    samtools view -bS - | \
    samtools sort -@ "$THREADS" -m 1G -T "$OUT_DIR/sort_tmp/sort" -o "$RAW_BAM" -

    samtools index -@ "$THREADS" "$RAW_BAM"
    rm -rf "$OUT_DIR/sort_tmp"
else
    echo "[STEP 2] Raw BAM already exists, skipping alignment."
fi

# ================= 2.5 PCR deduplication =================
DEDUP_BAM="$OUT_DIR/${PREFIX}.dedup.bam"
if [ ! -f "${DEDUP_BAM}.bai" ]; then
    echo "[STEP 2.5] Removing PCR duplicates with Sambamba..."
    mkdir -p "$OUT_DIR/sambamba_tmp"
    sambamba markdup -r -t 8 --tmpdir="$OUT_DIR/sambamba_tmp" "$RAW_BAM" "$DEDUP_BAM"
    rm -rf "$OUT_DIR/sambamba_tmp"
    samtools index -@ "$THREADS" "$DEDUP_BAM"
else
    echo "[STEP 2.5] Deduplicated BAM already exists, skipping."
fi

# ================= 3. Tn5 cut-site shift correction =================
SHIFTED_BAM="$OUT_DIR/${PREFIX}.dedup.shifted.bam"
if [ ! -f "${SHIFTED_BAM}.bai" ]; then
    echo "[STEP 3] Applying Tn5 shift correction (+4/-5 bp) with deepTools alignmentSieve..."
    alignmentSieve -b "$DEDUP_BAM" -o "${SHIFTED_BAM}.tmp" --ATACshift -p "$THREADS"
    samtools sort -@ "$THREADS" -o "$SHIFTED_BAM" "${SHIFTED_BAM}.tmp"
    samtools index -@ "$THREADS" "$SHIFTED_BAM"
    rm "${SHIFTED_BAM}.tmp"
else
    echo "[STEP 3] Shifted BAM already exists, skipping."
fi

# ================= 4. Peak calling with MACS2 =================
PEAK_FILE="$OUT_DIR/${PREFIX}_peaks.narrowPeak"
if [ ! -f "$PEAK_FILE" ]; then
    echo "[STEP 4] Calling peaks with MACS2..."
    macs2 callpeak \
        -t "$SHIFTED_BAM" \
        -f BAMPE \
        -g 1.0e10 \
        -n "$PREFIX" \
        --outdir "$OUT_DIR" \
        -q 0.05 \
        --keep-dup all
else
    echo "[STEP 4] Peak file already exists, skipping MACS2."
fi

echo "[SUCCESS] Phase 00 complete: alignment and peak calling finished."
