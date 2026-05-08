"""Parse all (rs.txt, ref.txt) sample pairs into joined parquets, fast and lean.

Memory-tight: scan_csv -> LazyFrame -> sink_parquet(engine='streaming').
The streaming engine batches parse + dedup + join + write, so per-worker
peak stays in the low GBs even for the largest 350 MB / 6.5 M-row samples.
ThreadPoolExecutor over samples; atomic rename writes; skips samples whose
parquet already exists.

Usage:
    python parse_all.py <rundir> [workers=8]
"""
from __future__ import annotations

import concurrent.futures
import os
import sys
from pathlib import Path

import polars as pl

INFO = [
    "count",
    "avg_mapping_quality",
    "avg_basequality",
    "avg_se_mapping_quality",
    "num_plus_strand",
    "num_minus_strand",
    "avg_pos_as_fraction",
    "avg_num_mismatches_as_fraction",
    "avg_sum_mismatch_qualities",
    "num_q2_containing_reads",
    "avg_distance_to_q2_start_in_q2_reads",
    "avg_clipped_length",
    "avg_distance_to_effective_3p_end",
]

EMPTY_SCHEMA = {
    "chrom": pl.Utf8,
    "pos": pl.Int64,
    "base": pl.Utf8,
    **{c: pl.Float64 for c in INFO},
}

# Samples whose upstream `<sample>.bamReadCount.txt` was generated against a
# different reference fasta than the GATK-bundle hg38 the rs run uses. At
# scattered positions on one or two chromosomes per sample, upstream's column-3
# reference base is `N` while rs reports the actual base — bam-readcount emits
# avg_sum_mismatch_qualities=0 for every base entry at those positions. See
# README "Caveats" #4. The IDs themselves are kept in `bench/excluded_samples.txt`
# (gitignored — cohort metadata is access-controlled). Missing file → empty set.
def _load_excluded() -> frozenset[str]:
    p = Path(__file__).parent / "excluded_samples.txt"
    if not p.exists():
        return frozenset()
    return frozenset(
        line.strip() for line in p.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    )


EXCLUDED_SAMPLES = _load_excluded()


def parse_brc_lazy(path: Path) -> pl.LazyFrame:
    # \x01 (SOH) is never present in bam-readcount text — using it as the CSV
    # separator gives one Utf8 cell per line, straight into Arrow buffers,
    # without going through a Python str + list of strs.
    return (
        pl.scan_csv(
            path,
            has_header=False,
            separator="\x01",
            new_columns=["line"],
            schema={"line": pl.Utf8},
            quote_char=None,
        )
        .with_columns(pl.col("line").str.split("\t").alias("parts"))
        .filter(pl.col("parts").list.len() >= 5)
        .select(
            pl.col("parts").list.get(0).alias("chrom"),
            pl.col("parts").list.get(1).cast(pl.Int64).alias("pos"),
            pl.col("parts").list.slice(4).alias("base_fields"),
        )
        .explode("base_fields")
        .with_columns(pl.col("base_fields").str.split(":").alias("vals"))
        .filter(pl.col("vals").list.len() >= len(INFO) + 1)
        .with_columns(pl.col("vals").list.get(0).alias("base"))
        .filter(pl.col("base") != "=")
        .with_columns(pl.col("vals").list.get(1).cast(pl.Int64).alias("_count"))
        .filter(pl.col("_count") > 0)
        .with_columns(*[
            pl.col("vals").list.get(i + 1).cast(pl.Float64).alias(name)
            for i, name in enumerate(INFO)
        ])
        .select(["chrom", "pos", "base"] + INFO)
    )


def parse_pair(sample_dir: Path, joined_dir: Path) -> str:
    sid = sample_dir.name
    if sid in EXCLUDED_SAMPLES:
        return f"exclude {sid} (ref-fasta-N — see README caveat 4)"
    out = joined_dir / f"{sid}.parquet"
    if out.exists():
        return f"skip {sid}"
    rs = sample_dir / "rs.txt"
    ref = sample_dir / "ref.txt"
    if not rs.exists():
        return f"no rs.txt {sid}"
    if not (ref.is_file() or ref.is_symlink()):
        return f"no ref.txt {sid}"
    if rs.stat().st_size == 0:
        return f"empty rs.txt {sid}"

    ref_lf = (
        parse_brc_lazy(ref)
        .unique(subset=["chrom", "pos", "base"], keep="first")
        .rename({c: f"ref_{c}" for c in INFO})
    )
    rs_lf = parse_brc_lazy(rs).rename({c: f"rs_{c}" for c in INFO})
    joined = ref_lf.join(rs_lf, on=["chrom", "pos", "base"], how="inner")

    tmp = joined_dir / f".{sid}.parquet.tmp.{os.getpid()}.{os.urandom(4).hex()}"
    joined.sink_parquet(tmp, compression="zstd", engine="streaming")
    tmp.rename(out)
    return f"{sid}: ok"


def main():
    rundir = Path(sys.argv[1])
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    raw = rundir / "raw"
    joined = rundir / "joined"
    joined.mkdir(exist_ok=True)

    sample_dirs = sorted(d for d in raw.iterdir() if d.is_dir())
    todo = [d for d in sample_dirs if not (joined / f"{d.name}.parquet").exists()]
    print(
        f"{len(sample_dirs)} samples total, "
        f"{len(sample_dirs) - len(todo)} already parsed, "
        f"{len(todo)} to do, workers={workers}",
        flush=True,
    )

    if not todo:
        print(f"all done. parquets: {sum(1 for _ in joined.glob('*.parquet'))}")
        return

    pl.Config.set_streaming_chunk_size(50_000)
    os.environ.setdefault("POLARS_MAX_THREADS", "2")
    done = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(parse_pair, d, joined): d for d in todo}
        for fut in concurrent.futures.as_completed(futs):
            done += 1
            try:
                msg = fut.result()
            except Exception as e:
                msg = f"FAIL {futs[fut].name}: {e}"
            if done % 25 == 0 or done == 1 or done == len(todo):
                print(f"  [{done}/{len(todo)}] {msg}", flush=True)

    print(f"done. parquets: {sum(1 for _ in joined.glob('*.parquet'))}")


if __name__ == "__main__":
    main()
