#!/usr/bin/env bash
# =============================================================================
# create_test_data.sh
#
# Creates subsampled test FASTQs for the MAGMA pipeline test profile.
# Downloads 3 M. tuberculosis H37Rv samples from the EXIT-RIF cohort
# (PRJNA1026351) via EBI FTP and subsamples to 200,000 read pairs each
# using fastp.
#
# Subsampling depth rationale:
#   Mtb genome: ~4.4 Mb | read length: 151 bp
#   200,000 pairs × 151 bp × 2 / 4,400,000 ≈ 13.7x median coverage
#   → safely clears MAGMA's cutoff_median_coverage = 10 threshold
#
# Output files (upload to Zenodo after running):
#   SRR26331590_sub_1.fastq.gz  SRR26331590_sub_2.fastq.gz
#   SRR26331595_sub_1.fastq.gz  SRR26331595_sub_2.fastq.gz
#   SRR26331599_sub_1.fastq.gz  SRR26331599_sub_2.fastq.gz
#
# After uploading to Zenodo:
#   1. Record the DOI and per-file download URLs.
#   2. Update assets/samplesheet_test.csv with the Zenodo URLs.
#   3. Update conf/test.config with the new Zenodo base URL if needed.
#
# Requirements:
#   - fastp  ≥ 0.23  (conda install -c bioconda fastp)
#   - wget or curl
#
# Usage:
#   bash scripts/create_test_data.sh [--outdir <dir>] [--threads <n>] [--pairs <n>]
#
# Options:
#   --outdir   Output directory  (default: ./test_data_subsampled)
#   --threads  fastp threads      (default: 4)
#   --pairs    Read pairs to keep (default: 200000)
# =============================================================================

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────────
OUTDIR="./test_data_subsampled"
THREADS=4
PAIRS=200000

# ── Argument parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)  OUTDIR="$2";  shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --pairs)   PAIRS="$2";   shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Samples ──────────────────────────────────────────────────────────────────
# Three M. tuberculosis H37Rv samples from EXIT-RIF (PRJNA1026351)
# These are clinical isolates used as MAGMA test data.
declare -A SAMPLES
SAMPLES["SRR26331590"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR263/090/SRR26331590"
SAMPLES["SRR26331595"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR263/095/SRR26331595"
SAMPLES["SRR26331599"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR263/099/SRR26331599"

# ── Preflight checks ─────────────────────────────────────────────────────────
if ! command -v fastp &>/dev/null; then
    echo "ERROR: fastp not found. Install with: conda install -c bioconda fastp"
    exit 1
fi

FASTP_VERSION=$(fastp --version 2>&1 | head -1)
echo "Using $FASTP_VERSION"
echo "Subsampling to ${PAIRS} read pairs per sample"
echo "Output directory: ${OUTDIR}"
echo ""

mkdir -p "${OUTDIR}"
TMPDIR="${OUTDIR}/.tmp"
mkdir -p "${TMPDIR}"

# ── Download + subsample ─────────────────────────────────────────────────────
for ACCESSION in "${!SAMPLES[@]}"; do
    BASE_URL="${SAMPLES[$ACCESSION]}"
    R1_URL="${BASE_URL}/${ACCESSION}_1.fastq.gz"
    R2_URL="${BASE_URL}/${ACCESSION}_2.fastq.gz"

    R1_RAW="${TMPDIR}/${ACCESSION}_1.fastq.gz"
    R2_RAW="${TMPDIR}/${ACCESSION}_2.fastq.gz"

    R1_SUB="${OUTDIR}/${ACCESSION}_sub_1.fastq.gz"
    R2_SUB="${OUTDIR}/${ACCESSION}_sub_2.fastq.gz"

    # Skip if already done
    if [[ -f "${R1_SUB}" && -f "${R2_SUB}" ]]; then
        echo "[${ACCESSION}] Subsampled files already present — skipping."
        continue
    fi

    # ── Download ────────────────────────────────────────────────────────────
    echo "[${ACCESSION}] Downloading R1..."
    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "${R1_RAW}" "${R1_URL}"
    else
        curl -# -L -o "${R1_RAW}" "${R1_URL}"
    fi

    echo "[${ACCESSION}] Downloading R2..."
    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "${R2_RAW}" "${R2_URL}"
    else
        curl -# -L -o "${R2_RAW}" "${R2_URL}"
    fi

    # ── Subsample with fastp ────────────────────────────────────────────────
    echo "[${ACCESSION}] Subsampling to ${PAIRS} read pairs..."
    fastp \
        --in1  "${R1_RAW}" \
        --in2  "${R2_RAW}" \
        --out1 "${R1_SUB}" \
        --out2 "${R2_SUB}" \
        --reads_to_process "${PAIRS}" \
        --disable_adapter_trimming \
        --disable_quality_filtering \
        --disable_length_filtering \
        --thread "${THREADS}" \
        --json "${OUTDIR}/${ACCESSION}_fastp.json" \
        --html "${OUTDIR}/${ACCESSION}_fastp.html" \
        2>&1 | grep -E "^(Read|Filtering|total|passed)" || true

    # ── Clean up raw download ───────────────────────────────────────────────
    rm -f "${R1_RAW}" "${R2_RAW}"

    R1_SIZE=$(du -h "${R1_SUB}" | cut -f1)
    R2_SIZE=$(du -h "${R2_SUB}" | cut -f1)
    echo "[${ACCESSION}] Done → R1: ${R1_SIZE}  R2: ${R2_SIZE}"
    echo ""
done

rmdir "${TMPDIR}" 2>/dev/null || true

# ── Summary ──────────────────────────────────────────────────────────────────
echo "========================================"
echo "Subsampled files ready in: ${OUTDIR}"
echo ""
ls -lh "${OUTDIR}"/*.fastq.gz 2>/dev/null
echo ""

TOTAL=$(du -sh "${OUTDIR}"/*.fastq.gz 2>/dev/null | tail -1 | cut -f1 || echo "?")
echo "Next steps:"
echo "  1. Upload the 6 .fastq.gz files to Zenodo"
echo "     (Create a new record, add all files, publish)"
echo ""
echo "  2. Record the per-file download URLs:"
echo "     https://zenodo.org/records/<RECORD_ID>/files/<FILENAME>?download=1"
echo ""
echo "  3. Update assets/samplesheet_test.csv with those URLs"
echo ""
echo "  4. Commit assets/samplesheet_test.csv (not the FASTQ files — they live on Zenodo)"
echo "========================================"
