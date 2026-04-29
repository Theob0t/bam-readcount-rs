"""Stream-parse a bam-readcount text file into a long-form parquet.

Memory-efficient pre-processing for the benchmark notebook: each sample's
ref+rs outputs get parsed once into a parquet file, joined on (chrom, pos,
base), and stored as bench/results/<runid>/joined/<sample>.parquet. The
notebook then scans those parquets via polars (lazy, low RAM).

Usage:
    python parse_brc.py <ref_brc> <rs_brc> <out_parquet>
"""
from __future__ import annotations

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
    """Vectorized polars parse of one bam-readcount text file.
    One row per (chrom, pos, base) where count > 0; '=' base skipped.
    """
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
    df = df.with_columns(metric_cols).select(["chrom", "pos", "base"] + INFO)
    return df


def main():
    ref_path, rs_path, out_path = (
        Path(sys.argv[1]),
        Path(sys.argv[2]),
        Path(sys.argv[3]),
    )
    ref_df = parse_brc(ref_path).unique(subset=["chrom", "pos", "base"], keep="first")
    ref_df = ref_df.rename({c: f"ref_{c}" for c in INFO})
    rs_df = parse_brc(rs_path).rename({c: f"rs_{c}" for c in INFO})
    joined = ref_df.join(rs_df, on=["chrom", "pos", "base"], how="inner")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    joined.write_parquet(out_path, compression="zstd")
    print(
        f"{ref_path.name}: ref={len(ref_df):,} rs={len(rs_df):,} joined={len(joined):,} -> {out_path}"
    )


if __name__ == "__main__":
    main()
