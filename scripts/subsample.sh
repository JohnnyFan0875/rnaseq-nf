#!/usr/bin/env bash
# =============================================================================
# subsample.sh - Subsample paired-end FASTQ files using seqtk
#
# Usage:
#   bash subsample.sh -i <input_dir> -o <output_dir> [options]
#
# Description:
#   Randomly subsample all paired-end FASTQ files (*_1.fastq.gz / *_2.fastq.gz)
#   in the input directory by a given fraction. The same seed is applied to
#   both R1 and R2 to ensure pairs remain in sync. Output files are gzip
#   compressed and written to the output directory with the same filenames.
#
# Options:
#   -i  <path>    Input directory containing *_1.fastq.gz and *_2.fastq.gz  [required]
#   -o  <path>    Output directory for subsampled FASTQ files               [required]
#   -f  <float>   Fraction of reads to keep, between 0 and 1 (default: 0.1)
#   -s  <int>     Random seed for reproducibility (default: 42)
#   -h            Show this help message and exit
#
# Requirements:
#   seqtk (https://github.com/lh3/seqtk)
#     Mac  : brew install seqtk
#     Linux: apt-get install seqtk
#     conda: conda install -c bioconda seqtk
#
# Examples:
#   # Subsample all samples at 10% (default)
#   bash subsample.sh -i data/GSE52778/raw -o data/GSE52778_sub/raw
#
#   # Subsample at 20% with a custom seed
#   bash subsample.sh -i data/GSE52778/raw -o data/GSE52778_sub/raw -f 0.2 -s 123
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT_DIR=""
OUTPUT_DIR=""
FRACTION="0.1"
SEED="42"

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    awk '/^# =/{found++} found==1{sub(/^# ?/,""); print} found==2{exit}' "$0"
    exit 0
}

# ── Argument parser ───────────────────────────────────────────────────────────
while getopts ":i:o:f:s:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        f) FRACTION="$OPTARG" ;;
        s) SEED="$OPTARG" ;;
        h) usage ;;
        :) echo "ERROR: Option -$OPTARG requires an argument." >&2; exit 1 ;;
        \?) echo "ERROR: Unknown option -$OPTARG" >&2; exit 1 ;;
    esac
done

# ── Validate required arguments ───────────────────────────────────────────────
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: -i and -o are required." >&2
    echo "Run 'bash $0 -h' for usage." >&2
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory not found: $INPUT_DIR" >&2
    exit 1
fi

if ! awk "BEGIN { exit !($FRACTION > 0 && $FRACTION < 1) }"; then
    echo "ERROR: -f must be between 0 and 1 (e.g. 0.1 for 10%)" >&2
    exit 1
fi

# ── Check seqtk ───────────────────────────────────────────────────────────────
if ! command -v seqtk &> /dev/null; then
    echo "ERROR: seqtk not found. Install it first:" >&2
    echo "  Mac  : brew install seqtk" >&2
    echo "  Linux: apt-get install seqtk" >&2
    echo "  conda: conda install -c bioconda seqtk" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "============================================"
echo " Input dir  : $INPUT_DIR"
echo " Output dir : $OUTPUT_DIR"
echo " Fraction   : $FRACTION"
echo " Seed       : $SEED"
echo " seqtk      : $(seqtk 2>&1 | head -2 | tail -1)"
echo "============================================"
echo ""

# ── Find R1 files ─────────────────────────────────────────────────────────────
R1_FILES=("$INPUT_DIR"/*_1.fastq.gz)

if [[ ${#R1_FILES[@]} -eq 0 || ! -f "${R1_FILES[0]}" ]]; then
    echo "ERROR: No *_1.fastq.gz files found in $INPUT_DIR" >&2
    exit 1
fi

echo "Found ${#R1_FILES[@]} sample(s)"
echo ""

# ── Subsample ─────────────────────────────────────────────────────────────────
for R1 in "${R1_FILES[@]}"; do

    BASENAME=$(basename "$R1" _1.fastq.gz)
    R2="${INPUT_DIR}/${BASENAME}_2.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "WARNING: R2 not found for $BASENAME, skipping..." >&2
        continue
    fi

    FRACTION_TAG=$(echo "$FRACTION" | sed 's/^\./0./')
    OUT_R1="${OUTPUT_DIR}/${BASENAME}_sub${FRACTION_TAG}_1.fastq.gz"
    OUT_R2="${OUTPUT_DIR}/${BASENAME}_sub${FRACTION_TAG}_2.fastq.gz"

    echo "Processing: $BASENAME"

    seqtk sample -s "$SEED" "$R1" "$FRACTION" | gzip > "$OUT_R1"
    seqtk sample -s "$SEED" "$R2" "$FRACTION" | gzip > "$OUT_R2"

    N=$(zcat "$OUT_R1" | awk 'NR % 4 == 1' | wc -l)
    echo "  -> ${N} reads retained"

done

echo ""
echo "Done! Subsampled files saved to: $OUTPUT_DIR"