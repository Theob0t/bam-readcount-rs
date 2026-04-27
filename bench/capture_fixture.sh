#!/usr/bin/env bash
# Capture a tiny region from a real BAM into a tests/data/fixtures/<name>/
# directory, and run upstream bam-readcount on it to produce the "expected.txt"
# oracle. Use this when chasing a benchmark discrepancy: extract the offending
# region, lock its expected output in, then iterate on the rs implementation
# until the integration test for that fixture passes.
#
# Usage:
#   bench/capture_fixture.sh <fixture_name> <bam> <region> <bed_file>
#   - fixture_name:  short identifier, e.g. D01_q2_reverse
#   - bam:           source BAM
#   - region:        samtools-style region string, e.g. chr1:924530-944953
#   - bed_file:      3-col BED to query (subset of positions in region)
set -euo pipefail

NAME=${1:?fixture name}
SRC_BAM=${2:?source bam}
REGION=${3:?region}
SRC_BED=${4:?bed file}

REPO=/gpfs/commons/home/tbotella/bam-readcount-rs
SAMTOOLS=/gpfs/commons/home/tbotella/miniconda3/envs/strega-hpc/bin/samtools
BRC=/gpfs/commons/home/tbotella/miniconda3/envs/strega-hpc/bin/bam-readcount
REF=/gpfs/commons/resources/GatkBundle/Homo_sapiens_assembly38.fasta

DIR=${REPO}/tests/data/fixtures/${NAME}
mkdir -p "$DIR"

cp "$SRC_BED" "$DIR/region.bed"
$SAMTOOLS view -b "$SRC_BAM" "$REGION" > "$DIR/region.bam" 2>/dev/null
$SAMTOOLS index "$DIR/region.bam"
$BRC -w1 -f "$REF" -l "$DIR/region.bed" "$DIR/region.bam" \
    > "$DIR/expected.txt" 2> "$DIR/brc.stderr"
echo "captured fixture '$NAME': $(wc -l < $DIR/expected.txt) expected lines"
echo "Add '$NAME' to FIXTURES in tests/integration.rs to enable the test."
