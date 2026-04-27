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

PY=/gpfs/commons/home/tbotella/miniconda3/envs/sobulk24/bin/papermill
$PY ${REPO}/bench/benchmark.ipynb ${RUNDIR}/benchmark.ipynb \
    -p rundir "${RUNDIR}" \
    --log-output 2>&1 | tee ${RUNDIR}/notebook.log

# Re-point latest symlink for README to embed
ln -sfn $(basename ${RUNDIR}) ${REPO}/bench/results/latest
echo "Done. Latest symlink -> ${REPO}/bench/results/latest"
