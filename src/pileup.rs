use crate::bed::Site;
use crate::metrics::{base_index, cigar_decode, PerBaseRecord, ReadScan, BASE_LABELS};
use anyhow::{Context, Result};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::faidx;

#[derive(Debug, Clone, Copy)]
pub struct PileupConfig {
    pub min_mapq: u8,
    pub min_baseq: u8,
    pub max_depth: u32,
}

pub struct PositionResult {
    pub chrom: String,
    pub pos: u32,
    pub ref_base: u8,
    pub depth: u32,
    pub per_base: [PerBaseRecord; 5],
}

/// Per-worker scratch buffers reused across chunks. Each worker owns one of
/// these; chunks are processed sequentially within a worker so reuse is safe.
pub struct ChunkScratch {
    /// Direct-indexed slot lookup: for an aligned 1-based position `p`,
    /// `slot_for_offset[p - want_min] = slot_idx` (or -1 if not wanted).
    /// Sized to (want_max - want_min + 1) per chunk, reused across chunks.
    slot_for_offset: Vec<i32>,
    /// Per-wanted-position accumulators (length = chunk.len()).
    slots: Vec<[PerBaseRecord; 5]>,
    /// Per-wanted-position depth (sum of all bases that passed filters).
    depths: Vec<u32>,
    /// One reusable Record into which we read from htslib (avoids per-record
    /// allocation).
    record: Record,
    /// Bulk reference cache for the current chunk's span.
    ref_cache: Vec<u8>,
}

impl ChunkScratch {
    pub fn new() -> Self {
        Self {
            slot_for_offset: Vec::new(),
            slots: Vec::new(),
            depths: Vec::new(),
            record: Record::new(),
            ref_cache: Vec::new(),
        }
    }
}

const SKIP_MASK: u16 = 0x4 | 0x100 | 0x200 | 0x400;

// CIGAR op codes (mirror metrics.rs constants — kept private there but BAM-spec
// stable, repeated here to keep the inner loop allocation-free without exposing
// internals).
const OP_MATCH: u8 = 0;
const OP_INS: u8 = 1;
const OP_DEL: u8 = 2;
const OP_REFSKIP: u8 = 3;
const OP_SOFT_CLIP: u8 = 4;
const OP_EQUAL: u8 = 7;
const OP_DIFF: u8 = 8;

/// Process one chunk (consecutive sites within one chrom) using a read-driven
/// walk: iterate records in the fetch range, compute ReadScan once per read,
/// then walk CIGAR M/=/X positions and emit observations into per-position
/// accumulator slots. Replaces the previous pileup-driven loop which
/// recomputed ReadScan once per (read × wanted-position-the-read-covers) and
/// paid full pileup-column overhead.
pub fn process_chunk(
    reader: &mut IndexedReader,
    fasta: &faidx::Reader,
    sites: &[Site],
    config: PileupConfig,
    scratch: &mut ChunkScratch,
) -> Result<Vec<PositionResult>> {
    if sites.is_empty() {
        return Ok(Vec::new());
    }

    let header = reader.header().clone();
    let chrom = sites[0].chrom.clone();
    let tid = match header.tid(chrom.as_bytes()) {
        Some(t) => t,
        None => return Ok(Vec::new()),
    };

    let want_min = sites.iter().map(|s| s.pos).min().unwrap();
    let want_max = sites.iter().map(|s| s.pos).max().unwrap();
    let span = (want_max - want_min + 1) as usize;

    scratch.slot_for_offset.clear();
    scratch.slot_for_offset.resize(span, -1);
    scratch.slots.clear();
    scratch.slots.resize(sites.len(), Default::default());
    scratch.depths.clear();
    scratch.depths.resize(sites.len(), 0);

    for (idx, site) in sites.iter().enumerate() {
        let off = (site.pos - want_min) as usize;
        // If the BED has duplicate positions, last one wins; benign because
        // sites with the same pos produce the same output.
        scratch.slot_for_offset[off] = idx as i32;
    }

    // bam-readcount uses 1-based positions on input (BED-like file); BAM is
    // 0-based half-open. Fetch [want_min - 1, want_max) in 0-based half-open
    // so the iterator returns every read overlapping any wanted position.
    let fetch_start = (want_min.saturating_sub(1)) as i64;
    let fetch_end = want_max as i64;
    reader
        .fetch((tid, fetch_start, fetch_end))
        .with_context(|| format!("fetching {}:{}-{}", chrom, fetch_start, fetch_end))?;

    // One bulk reference fetch for the chunk's span (vs N tiny fetches in the
    // pileup-driven version). faidx::fetch_seq is end-inclusive.
    let raw_ref = fasta
        .fetch_seq(&chrom, fetch_start as usize, (fetch_end - 1) as usize)
        .with_context(|| format!("fetching ref {}:{}-{}", chrom, fetch_start, fetch_end))?;
    scratch.ref_cache.clear();
    scratch.ref_cache.reserve(raw_ref.len());
    scratch
        .ref_cache
        .extend(raw_ref.iter().map(|b| b.to_ascii_uppercase()));

    // Per-record loop. `reader.read(&mut record)` reuses the record buffer.
    while let Some(result) = reader.read(&mut scratch.record) {
        result?;
        process_record(
            &scratch.record,
            want_min,
            want_max,
            &scratch.slot_for_offset,
            &mut scratch.slots,
            &mut scratch.depths,
            &config,
        );
    }

    // Build results in BED order. `sites` is already sorted by caller.
    let mut results: Vec<PositionResult> = Vec::with_capacity(sites.len());
    for (idx, site) in sites.iter().enumerate() {
        let ref_off = (site.pos as i64 - 1 - fetch_start) as usize;
        let ref_base = scratch
            .ref_cache
            .get(ref_off)
            .copied()
            .unwrap_or(b'N');
        let per_base = std::mem::take(&mut scratch.slots[idx]);
        let depth = scratch.depths[idx];
        results.push(PositionResult {
            chrom: chrom.clone(),
            pos: site.pos,
            ref_base,
            depth,
            per_base,
        });
    }

    Ok(results)
}

#[inline]
fn process_record(
    record: &Record,
    want_min: u32,
    want_max: u32,
    slot_for_offset: &[i32],
    slots: &mut [[PerBaseRecord; 5]],
    depths: &mut [u32],
    config: &PileupConfig,
) {
    // Standard mpileup filter mask (matches bam-readcount v1.0.1 which uses
    // BAM_DEF_MASK = FUNMAP|FSECONDARY|FQCFAIL|FDUP). Supplementary
    // alignments (0x800) are NOT excluded by upstream.
    if (record.flags() & SKIP_MASK) != 0 {
        return;
    }
    if record.mapq() < config.min_mapq {
        return;
    }
    let pos_0 = record.pos();
    if pos_0 < 0 {
        return;
    }
    let read_start_0 = pos_0 as u32;

    // Bounding-box pre-check: walk raw_cigar once to compute the read's
    // reference end, skip the read entirely if it doesn't overlap the
    // chunk's wanted span. For STREGA's clustered exonic BEDs this rarely
    // fires, but it's nearly free (one cigar pass).
    let raw = record.raw_cigar();
    let mut ref_len: u32 = 0;
    for &packed in raw {
        let (op, len) = cigar_decode(packed);
        match op {
            OP_MATCH | OP_EQUAL | OP_DIFF | OP_DEL | OP_REFSKIP => ref_len += len,
            _ => {}
        }
    }
    let read_end_1 = read_start_0 + ref_len; // 1-based last covered position
    let read_start_1 = read_start_0 + 1;
    if read_end_1 < want_min || read_start_1 > want_max {
        return;
    }

    // ReadScan is computed once per read (vs once per pileup hit before).
    let scan = ReadScan::from_record(record);

    let qual = record.qual();
    let seq = record.seq();
    let mut read_idx: u32 = 0;
    let mut ref_pos_0: u32 = read_start_0;

    for &packed in raw {
        let (op, len) = cigar_decode(packed);
        match op {
            OP_MATCH | OP_EQUAL | OP_DIFF => {
                for _ in 0..len {
                    let one_based = ref_pos_0 + 1;
                    if one_based >= want_min && one_based <= want_max {
                        let off = (one_based - want_min) as usize;
                        let slot_idx = slot_for_offset[off];
                        if slot_idx >= 0 {
                            let slot_idx = slot_idx as usize;
                            if depths[slot_idx] < config.max_depth {
                                let base_qual = qual[read_idx as usize];
                                if base_qual >= config.min_baseq {
                                    let read_base = seq[read_idx as usize];
                                    let obs = scan.observation_at(read_idx, base_qual);
                                    let bi = base_index(read_base);
                                    slots[slot_idx][bi].push(&obs);
                                    depths[slot_idx] += 1;
                                }
                            }
                        }
                    }
                    read_idx += 1;
                    ref_pos_0 += 1;
                }
            }
            OP_INS | OP_SOFT_CLIP => {
                read_idx += len;
            }
            OP_DEL | OP_REFSKIP => {
                ref_pos_0 += len;
            }
            // HardClip / Pad: no advance.
            _ => {}
        }
    }
}

pub fn format_position(result: &PositionResult, out: &mut String) {
    use std::fmt::Write;
    write!(
        out,
        "{}\t{}\t{}\t{}",
        result.chrom,
        result.pos,
        result.ref_base as char,
        result.depth
    )
    .unwrap();
    out.push('\t');
    // Always-zero "=" base record (matches upstream output framing).
    out.push_str("=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00");
    for (i, label) in BASE_LABELS.iter().enumerate() {
        out.push('\t');
        result.per_base[i].write_with_label(label, out);
    }
    out.push('\n');
}
