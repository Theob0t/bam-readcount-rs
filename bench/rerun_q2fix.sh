#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --account=dllab
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --job-name=brc_rerun
#SBATCH --output=/gpfs/commons/home/tbotella/bam-readcount-rs/bench/logs/%x_%j.log

# Re-run the 2000-sample benchmark on the q2-fixed rs.txt outputs:
#   1. cleanup: drop stale pre-q2-fix run, promote q2fix into canonical name,
#      wipe partial joined parquets, drop superseded 200-sample runs.
#   2. parse: parse_all.py over the 1955 raw dirs (skips the 17 in
#      EXCLUDED_SAMPLES via README "Caveats" #4).
#   3. notebook: papermill benchmark.ipynb to regenerate SUMMARY.md,
#      per_feature_corr.tsv, and the plots.
#
# Submit with:
#   mkdir -p /gpfs/commons/home/tbotella/bam-readcount-rs/bench/logs
#   sbatch /gpfs/commons/home/tbotella/bam-readcount-rs/bench/rerun_q2fix.sh

set -euo pipefail

REPO=/gpfs/commons/home/tbotella/bam-readcount-rs
RESULTS=$REPO/bench/results
PY=/gpfs/commons/home/tbotella/miniconda3/envs/sobulk24/bin/python
PMILL=/gpfs/commons/home/tbotella/miniconda3/envs/sobulk24/bin/papermill

echo "[$(date)] === cleanup ==="

# 1. Drop stale pre-q2-fix 2000samples/ (raw/ is from Apr 28, before cdfe28f).
rm -rf "$RESULTS/2000samples"

# 2. Promote q2fix into the canonical 2000samples name so README paths still work.
if [ -d "$RESULTS/2000samples_q2fix" ]; then
    mv "$RESULTS/2000samples_q2fix" "$RESULTS/2000samples"
fi

# 3. Wipe partial joined parquets + stale logs from the OOM-killed parse attempt.
rm -rf "$RESULTS/2000samples/joined"
rm -f  "$RESULTS/2000samples/parse.log"
rm -f  "$RESULTS/2000samples/logs/agg.log"

# 4. Repoint latest pointers.
ln -sfn 2000samples "$RESULTS/latest"
echo 2000samples > "$RESULTS/.latest_runid"

# 5. Drop superseded 200-sample runs (no tracked files in either).
rm -rf "$RESULTS/200samples_v1"
rm -rf "$RESULTS/200samples_v2"

echo "[$(date)] cleanup done. results dir state:"
ls -lh "$RESULTS/"
echo "raw subdirs: $(ls "$RESULTS/2000samples/raw" | wc -l)"

echo "[$(date)] === parse_all.py (8 workers x 2 polars threads = 16 total) ==="
$PY $REPO/bench/parse_all.py "$RESULTS/2000samples" 8 \
    2>&1 | tee "$RESULTS/2000samples/parse.log"

echo "[$(date)] joined parquets: $(ls "$RESULTS/2000samples/joined" | wc -l)"

echo "[$(date)] === papermill benchmark.ipynb ==="
$PMILL $REPO/bench/benchmark.ipynb "$RESULTS/2000samples/benchmark.ipynb" \
    -p rundir "$RESULTS/2000samples" \
    --log-output \
    2>&1 | tee "$RESULTS/2000samples/notebook.log"

echo "[$(date)] === done ==="
ls -lh "$RESULTS/2000samples/"
