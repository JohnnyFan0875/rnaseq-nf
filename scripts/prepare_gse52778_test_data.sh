#!/usr/bin/env bash
# =============================================================================
# prepare_gse52778_test_data.sh
#
# Download a small balanced subset of GSE52778 and subsample it for pipeline
# testing. The script keeps 2 control samples and 2 treatment samples, then
# writes a ready-to-use metadata CSV for the generated FASTQ files.
#
# Usage:
#   bash scripts/prepare_gse52778_test_data.sh -o data/GSE52778_test -f 0.05
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
#   curl
#   wget
#   seqtk
# =============================================================================

set -euo pipefail

OUTDIR=""
FRACTION=""
SEED="42"
ENA_STUDY="SRP033351"
ENA_API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ENA_STUDY}&result=read_run&fields=run_accession,fastq_ftp&format=tsv"

usage() {
    sed -n '2,23p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
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

for cmd in curl wget seqtk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: $cmd" >&2
        exit 1
    fi
done

RAW_DIR="${OUTDIR}/raw"
META_DIR="${OUTDIR}/metadata"
TMP_DIR="${OUTDIR}/tmp_download"
METADATA_FILE="${META_DIR}/sample_metadata.csv"
FRACTION_TAG="$(printf '%s' "$FRACTION" | sed 's/^\./0./')"

mkdir -p "$RAW_DIR" "$META_DIR" "$TMP_DIR"

SAMPLES=(
    "SRR1039508 control 1"
    "SRR1039509 treatment 1"
    "SRR1039512 control 2"
    "SRR1039513 treatment 2"
)

declare -A SAMPLE_GROUP=()
declare -A SAMPLE_REPLICATE=()
declare -A SAMPLE_R1_URL=()
declare -A SAMPLE_R2_URL=()

for entry in "${SAMPLES[@]}"; do
    read -r sample_id group replicate <<<"$entry"
    SAMPLE_GROUP["$sample_id"]="$group"
    SAMPLE_REPLICATE["$sample_id"]="$replicate"
done

echo "Fetching FASTQ locations from ENA..."
while IFS=$'\t' read -r run_accession fastq_ftp; do
    [[ -z "$run_accession" || "$run_accession" == "run_accession" ]] && continue
    [[ -z "${SAMPLE_GROUP[$run_accession]:-}" ]] && continue

    mapfile -t paths < <(printf '%s\n' "$fastq_ftp" | tr ';' '\n' | grep -E '_[12]\.fastq\.gz$')

    if [[ "${#paths[@]}" -ne 2 ]]; then
        echo "ERROR: Expected 2 FASTQ paths for ${run_accession}, found ${#paths[@]}." >&2
        exit 1
    fi

    for path in "${paths[@]}"; do
        if [[ "$path" == *_1.fastq.gz ]]; then
            SAMPLE_R1_URL["$run_accession"]="ftp://${path}"
        elif [[ "$path" == *_2.fastq.gz ]]; then
            SAMPLE_R2_URL["$run_accession"]="ftp://${path}"
        fi
    done
done < <(curl -fsSL "$ENA_API_URL")

printf 'sample_id,fastq1,fastq2,group,replicate\n' > "$METADATA_FILE"

for entry in "${SAMPLES[@]}"; do
    read -r sample_id group replicate <<<"$entry"

    r1_url="${SAMPLE_R1_URL[$sample_id]:-}"
    r2_url="${SAMPLE_R2_URL[$sample_id]:-}"

    if [[ -z "$r1_url" || -z "$r2_url" ]]; then
        echo "ERROR: Missing FASTQ URL(s) for ${sample_id}." >&2
        exit 1
    fi

    tmp_r1="${TMP_DIR}/${sample_id}_1.fastq.gz"
    tmp_r2="${TMP_DIR}/${sample_id}_2.fastq.gz"
    out_sample_id="${sample_id}_sub${FRACTION_TAG}"
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
