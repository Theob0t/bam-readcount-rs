#!/usr/bin/env bash
# Aggregate per-sample metrics + run the analysis notebook.
set -euo pipefail

RUNDIR=${1:?results dir}
REPO=/gpfs/commons/home/tbotella/bam-readcount-rs

HEADER_DONE=0
> ${RUNDIR}/per_sample_metrics.tsv
for f in ${RUNDIR}/raw/*/metrics.tsv; do
    [ -f "$f" ] || continue
    if [ "$HEADER_DONE" = "0" ]; then
        head -1 "$f" >> ${RUNDIR}/per_sample_metrics.tsv
        HEADER_DONE=1
    fi
    tail -n +2 "$f" >> ${RUNDIR}/per_sample_metrics.tsv
done
N=$(($(wc -l < ${RUNDIR}/per_sample_metrics.tsv) - 1))
echo "Aggregated $N rows -> ${RUNDIR}/per_sample_metrics.tsv"

PY=/gpfs/commons/home/tbotella/miniconda3/envs/sobulk24/bin/python

# Pre-parse each sample's (ref, rs) into a joined parquet. This is the heavy
# step (text parsing + join), done once per sample, in parallel, so the
# notebook only needs to scan parquets.
mkdir -p ${RUNDIR}/joined
for d in ${RUNDIR}/raw/*; do
    sid=$(basename "$d")
    out=${RUNDIR}/joined/${sid}.parquet
    [ -f "$out" ] && continue
    [ -s "$d/rs.txt" ] || continue
    [ -L "$d/ref.txt" ] || [ -f "$d/ref.txt" ] || continue
    echo "  parsing $sid..."
    $PY ${REPO}/bench/parse_brc.py "$d/ref.txt" "$d/rs.txt" "$out" 2>&1 \
        | tee -a ${RUNDIR}/parse.log &
    if (( $(jobs -r | wc -l) >= ${PARALLEL:-2} )); then
        wait -n
    fi
done
wait
echo "all parses done. parquets: $(ls ${RUNDIR}/joined/ | wc -l)"

PMILL=/gpfs/commons/home/tbotella/miniconda3/envs/sobulk24/bin/papermill
$PMILL ${REPO}/bench/benchmark.ipynb ${RUNDIR}/benchmark.ipynb \
    -p rundir "${RUNDIR}" \
    --log-output 2>&1 | tee ${RUNDIR}/notebook.log

ln -sfn $(basename ${RUNDIR}) ${REPO}/bench/results/latest
echo "Done. Latest symlink -> ${REPO}/bench/results/latest"
