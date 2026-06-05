#!/usr/bin/env bash
# =============================================================================
# prepare_test_data.sh
#
# Download a small balanced subset of GSE52778 and subsample it for pipeline
# testing. The script keeps 2 control samples and 2 treatment samples, then
# writes a ready-to-use metadata CSV for the generated FASTQ files.
#
# Usage:
#   bash scripts/prepare_test_data.sh -o data/GSE52778_test -f 0.05
#
# Options:
#   -o  <path>   Output dataset directory                               [required]
#   -f  <float>  Fraction of reads to keep, between 0 and 1             [required]
#   -h           Show this help message and exit
#
# Output structure:
#   <outdir>/
#   ├── metadata/sample_metadata.csv
#   └── raw/*.fastq.gz
#
# Requirements:
#   wget
#   seqtk
# =============================================================================

set -euo pipefail

OUTDIR=""
FRACTION=""
SEED="42"

usage() {
    sed -n '2,23p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
}

# ---------------------------------------------------------------------------
# Compute EBI FTP base URL from a SRR accession.
# Rule: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/<first-6-chars>/00<last-digit>/<srr>
# ---------------------------------------------------------------------------
ena_ftp_base() {
    local srr="$1"
    local prefix="${srr:0:6}"
    local suffix="${srr: -1}"
    echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/00${suffix}/${srr}"
}

while getopts ":o:f:h" opt; do
    case "$opt" in
        o) OUTDIR="$OPTARG" ;;
        f) FRACTION="$OPTARG" ;;
        h) usage ;;
        :) echo "ERROR: Option -$OPTARG requires an argument." >&2; exit 1 ;;
        \?) echo "ERROR: Unknown option -$OPTARG" >&2; exit 1 ;;
    esac
done

if [[ -z "$OUTDIR" || -z "$FRACTION" ]]; then
    echo "ERROR: -o and -f are required." >&2
    echo "Run 'bash $0 -h' for usage." >&2
    exit 1
fi

if ! awk "BEGIN { exit !($FRACTION > 0 && $FRACTION < 1) }"; then
    echo "ERROR: -f must be between 0 and 1 (for example 0.05)." >&2
    exit 1
fi

for cmd in wget seqtk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: $cmd" >&2
        exit 1
    fi
done

RAW_DIR="${OUTDIR}/raw"
META_DIR="${OUTDIR}/metadata"
TMP_DIR="${OUTDIR}/tmp_download"
METADATA_FILE="${META_DIR}/sample_metadata.csv"

mkdir -p "$RAW_DIR" "$META_DIR" "$TMP_DIR"

# ---------------------------------------------------------------------------
# Sample definitions: "<SRR> <group> <replicate>"
#   SRR1039508  GSM1275862  N61311_untreated
#   SRR1039509  GSM1275863  N61311_Dex
#   SRR1039512  GSM1275866  N052611_untreated
#   SRR1039513  GSM1275867  N052611_Dex
# ---------------------------------------------------------------------------
SAMPLES=(
    "SRR1039508 untreated 1"
    "SRR1039509 dex       1"
    "SRR1039512 untreated 2"
    "SRR1039513 dex       2"
)

printf 'sample_id,fastq1,fastq2,group,replicate\n' > "$METADATA_FILE"

for entry in "${SAMPLES[@]}"; do
    read -r sample_id group replicate <<<"$entry"

    base_url="$(ena_ftp_base "$sample_id")"
    r1_url="${base_url}/${sample_id}_1.fastq.gz"
    r2_url="${base_url}/${sample_id}_2.fastq.gz"

    tmp_r1="${TMP_DIR}/${sample_id}_1.fastq.gz"
    tmp_r2="${TMP_DIR}/${sample_id}_2.fastq.gz"
    out_sample_id="${sample_id}_sub${FRACTION}"
    out_r1="${RAW_DIR}/${out_sample_id}_1.fastq.gz"
    out_r2="${RAW_DIR}/${out_sample_id}_2.fastq.gz"

    echo "Downloading ${sample_id}..."
    wget -O "$tmp_r1" "$r1_url"
    wget -O "$tmp_r2" "$r2_url"

    echo "Subsampling ${sample_id} at fraction ${FRACTION}..."
    seqtk sample -s "$SEED" "$tmp_r1" "$FRACTION" | gzip > "$out_r1"
    seqtk sample -s "$SEED" "$tmp_r2" "$FRACTION" | gzip > "$out_r2"

    rm -f "$tmp_r1" "$tmp_r2"

    printf '%s,%s,%s,%s,%s\n' \
        "$out_sample_id" \
        "$out_r1" \
        "$out_r2" \
        "$group" \
        "$replicate" >> "$METADATA_FILE"
done

rmdir "$TMP_DIR" 2>/dev/null || true

echo ""
echo "Done."
echo "Dataset directory : $OUTDIR"
echo "Metadata file     : $METADATA_FILE"
echo "Subsample fraction: $FRACTION"