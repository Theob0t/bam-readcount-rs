#!/usr/bin/env bash
# Run the benchmark notebook locally (no Slurm). Use after the array +
# parse_all have finished and joined/*.parquet is complete. The notebook
# does a full streaming scan over ~2.1B rows; on the login node polars
# streaming peak stays in low-GB.
set -euo pipefail

RUNID=${1:-2000samples_l2l3}
REPO=/gpfs/commons/home/tbotella/bam-readcount-rs
RUNDIR=${REPO}/bench/results/${RUNID}
PMILL=/gpfs/commons/home/tbotella/miniconda3/envs/sobulk24/bin/papermill

mkdir -p "${RUNDIR}/logs"

# Refuse to run while parse_all.py is still finalizing a sample — that race
# (partial <sid>.parquet visible to the scan) is what broke the previous run.
TMPS=$(find "${RUNDIR}/joined" -maxdepth 1 -name '.*.parquet.tmp.*' 2>/dev/null | wc -l)
if [ "${TMPS}" -gt 0 ]; then
    echo "ERROR: ${TMPS} in-flight tmp parquet(s) in ${RUNDIR}/joined — wait for parse_all to finish." >&2
    exit 1
fi

echo "[$(date)] papermill -> ${RUNDIR}/notebook.log"
${PMILL} ${REPO}/bench/benchmark.ipynb ${RUNDIR}/benchmark.ipynb \
    -p rundir "${RUNDIR}" \
    --log-output \
    2>&1 | tee "${RUNDIR}/notebook.log"

echo "[$(date)] done. Output: ${RUNDIR}/SUMMARY.md"
