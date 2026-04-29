"""Parse all (rs.txt, ref.txt) sample pairs into joined parquets, fast.

One Python process, ThreadPoolExecutor, vectorized polars parser, atomic
rename writes. Skips samples whose parquet already exists. Safe to run
concurrently with bench/aggregate.sh — the rename is atomic and the
parquets are byte-identical anyway.

Usage:
    python parse_all.py <rundir> [workers=16]
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


def parse_brc(path: Path) -> pl.DataFrame:
    text = Path(path).read_text()
    if not text.strip():
        return pl.DataFrame(schema=EMPTY_SCHEMA)
    lines = text.rstrip("\n").split("\n")
    df = pl.DataFrame({"line": lines})
    df = df.with_columns(pl.col("line").str.split("\t").alias("parts"))
    df = df.filter(pl.col("parts").list.len() >= 5)
    df = df.select(
        pl.col("parts").list.get(0).alias("chrom"),
        pl.col("parts").list.get(1).cast(pl.Int64).alias("pos"),
        pl.col("parts").list.slice(4).alias("base_fields"),
    ).explode("base_fields")
    df = df.with_columns(pl.col("base_fields").str.split(":").alias("vals"))
    df = df.filter(pl.col("vals").list.len() >= len(INFO) + 1)
    df = df.with_columns(pl.col("vals").list.get(0).alias("base"))
    df = df.filter(pl.col("base") != "=")
    df = df.with_columns(pl.col("vals").list.get(1).cast(pl.Int64).alias("_count"))
    df = df.filter(pl.col("_count") > 0)
    metric_cols = [
        pl.col("vals").list.get(i + 1).cast(pl.Float64).alias(name)
        for i, name in enumerate(INFO)
    ]
    return df.with_columns(metric_cols).select(["chrom", "pos", "base"] + INFO)


def parse_pair(sample_dir: Path, joined_dir: Path) -> str:
    sid = sample_dir.name
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
    ref_df = parse_brc(ref).unique(subset=["chrom", "pos", "base"], keep="first")
    ref_df = ref_df.rename({c: f"ref_{c}" for c in INFO})
    rs_df = parse_brc(rs).rename({c: f"rs_{c}" for c in INFO})
    joined = ref_df.join(rs_df, on=["chrom", "pos", "base"], how="inner")
    tmp = joined_dir / f".{sid}.parquet.tmp.{os.getpid()}.{os.urandom(4).hex()}"
    joined.write_parquet(tmp, compression="zstd")
    tmp.rename(out)
    return f"{sid}: ref={len(ref_df):,} rs={len(rs_df):,} joined={len(joined):,}"


def main():
    rundir = Path(sys.argv[1])
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 16
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
