#!/usr/bin/env bash
# Script        : run.sh
# Description   : Download the PanG2P test data, build the input CSVs and run pang2p.
# Date          : 2026-11-23
# Version       : 0.1.0
# Contact       : c.s.sivasubramani@gmail.com
#
# Usage:
#   ./run.sh                      # run pang2p with default settings
#   ./run.sh --help               # show pang2p help
#   THREADS=8 MEM=16g ./run.sh    # customize resources
#   DRY_RUN=1 ./run.sh            # only show what would run


set -euo pipefail

URL="https://github.com/tamu-plant-genomics/PanG2P/releases/download/v0.1.0-data/test_data.tar.gz"
WORKDIR="${WORKDIR:-$PWD/pang2p_test}"
THREADS="${THREADS:-4}"
MEM="${MEM:-8g}"
OUTDIR="${OUTDIR:-$WORKDIR/results}"

command -v pang2p >/dev/null 2>&1 || {
    echo "[run.sh] pang2p not found on PATH. Please install pang2p and ensure it is available in your PATH." >&2
    exit 1
}

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# ---- 1. download and unpack the test data (skipped if already present)
if [[ ! -f test_data.tar.gz ]]; then
    echo "[run.sh] downloading test data ..."
    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 3 -o test_data.tar.gz "$URL"
    elif command -v wget >/dev/null 2>&1; then
        wget -O test_data.tar.gz "$URL"
    else
        echo "[run.sh] need curl or wget" >&2
        exit 1
    fi
fi

if [[ ! -d test_data ]]; then
    echo "[run.sh] extracting ..."
    tar -xzf test_data.tar.gz
fi

GENOME_DIR=$(find "$WORKDIR" -type d -name genomes | head -n 1)
SAMPLE_DIR=$(find "$WORKDIR" -type d -name samples | head -n 1)
if [[ -z "$GENOME_DIR" || -z "$SAMPLE_DIR" ]]; then
    echo "[run.sh] could not find genomes/ and samples/ after extraction" >&2
    exit 1
fi

# ---- 2. generate the CSVs (absolute paths, ids taken from the file names)
echo "id,path" > genomes.csv
for f in "$GENOME_DIR"/*.fa "$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fna; do
    [[ -e "$f" ]] || continue
    id=$(basename "$f"); id=${id%.*}
    echo "$id,$(realpath "$f")" >> genomes.csv
done

echo "sample_id,path" > samples.csv
for f in "$SAMPLE_DIR"/*.fq.gz "$SAMPLE_DIR"/*.fastq.gz; do
    [[ -e "$f" ]] || continue
    id=$(basename "$f"); id=${id%.fq.gz}; id=${id%.fastq.gz}
    echo "$id,$(realpath "$f")" >> samples.csv
done

n_g=$(( $(wc -l < genomes.csv) - 1 ))
n_s=$(( $(wc -l < samples.csv) - 1 ))
echo "[run.sh] $n_g genomes, $n_s samples"
(( n_g > 0 && n_s > 0 )) || { echo "[run.sh] no inputs found" >&2; exit 1; }

# ---- 3. run the pipeline
EXTRA=()
[[ "${DRY_RUN:-0}" == "1" ]] && EXTRA+=(--dry-run)

pang2p \
    -g genomes.csv \
    -s samples.csv \
    -t "$THREADS" \
    -m "$MEM" \
    -o "$OUTDIR" \
    "${EXTRA[@]}" \
    "$@"

echo "[run.sh] done. Results: $OUTDIR"
echo "[run.sh] genotype table: $OUTDIR/work/04.genotyping/final.genotype.tsv"
