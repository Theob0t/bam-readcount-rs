use crate::bed::Site;
use crate::metrics::{base_index, PerBaseRecord, ReadScan, BASE_LABELS};
use anyhow::{Context, Result};
use rust_htslib::bam::{self, IndexedReader, Read};
use rust_htslib::faidx;
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

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

/// Process one chunk (consecutive sites within one chrom). Opens its own
/// IndexedReader so this can be called from multiple worker threads.
pub fn process_chunk(
    bam_path: &Path,
    fasta_path: &Path,
    sites: &[Site],
    config: PileupConfig,
) -> Result<Vec<PositionResult>> {
    if sites.is_empty() {
        return Ok(Vec::new());
    }

    let mut reader = IndexedReader::from_path(bam_path)
        .with_context(|| format!("opening BAM {}", bam_path.display()))?;
    let header = reader.header().clone();
    let chrom = sites[0].chrom.clone();
    let tid = match header.tid(chrom.as_bytes()) {
        Some(t) => t,
        None => return Ok(Vec::new()),
    };

    let min_pos = sites.iter().map(|s| s.pos).min().unwrap();
    let max_pos = sites.iter().map(|s| s.pos).max().unwrap();
    // bam-readcount uses 1-based positions on input (BED-like file); BAM is
    // 0-based half-open. Fetch [min_pos-1, max_pos) so the pileup covers all
    // queried sites. Use a slight overscan to capture reads that span the
    // first position.
    let fetch_start = (min_pos.saturating_sub(1)) as i64;
    let fetch_end = max_pos as i64;
    reader
        .fetch((tid, fetch_start, fetch_end))
        .with_context(|| format!("fetching {}:{}-{}", chrom, fetch_start, fetch_end))?;

    let fasta = faidx::Reader::from_path(fasta_path)
        .with_context(|| format!("opening FASTA {}", fasta_path.display()))?;

    // Build a lookup of which 1-based positions on this chrom we want to emit.
    let want: HashMap<u32, ()> = sites.iter().map(|s| (s.pos, ())).collect();
    let want_min = min_pos;
    let want_max = max_pos;

    // Per-read scan cache (keyed by qname-hash) so the MD/Q2/SA work is done
    // once per read even if the read spans multiple queried positions.
    let mut scan_cache: HashMap<u64, Arc<ReadScan>> = HashMap::new();

    let mut results: Vec<PositionResult> = Vec::with_capacity(sites.len());

    let mut pileup = reader.pileup();
    pileup.set_max_depth(config.max_depth);

    while let Some(p) = pileup.next() {
        let column = p?;
        let zero_based = column.pos();
        let one_based = zero_based + 1;
        if one_based < want_min || one_based > want_max {
            continue;
        }
        if !want.contains_key(&one_based) {
            continue;
        }

        let ref_base = fetch_ref_base(&fasta, &chrom, zero_based)?;
        let mut per_base = [
            PerBaseRecord::new(),
            PerBaseRecord::new(),
            PerBaseRecord::new(),
            PerBaseRecord::new(),
            PerBaseRecord::new(),
        ];
        let mut depth: u32 = 0;

        for alignment in column.alignments() {
            if alignment.is_del() || alignment.is_refskip() {
                continue;
            }
            let record = alignment.record();
            // Standard mpileup filter mask (matches bam-readcount v1.0.1
            // which uses BAM_DEF_MASK = FUNMAP|FSECONDARY|FQCFAIL|FDUP).
            // Supplementary alignments (0x800) are NOT excluded by upstream.
            const SKIP_MASK: u16 = 0x4 | 0x100 | 0x200 | 0x400;
            if (record.flags() & SKIP_MASK) != 0 {
                continue;
            }
            if record.mapq() < config.min_mapq {
                continue;
            }
            let qpos = match alignment.qpos() {
                Some(q) => q as u32,
                None => continue,
            };
            let seq = record.seq();
            let read_base = seq[qpos as usize];
            let qual = record.qual();
            let base_qual = qual[qpos as usize];
            if base_qual < config.min_baseq {
                continue;
            }

            // Get-or-build the per-read scan.
            let key = qname_hash(record.qname(), record.flags());
            let scan = scan_cache
                .entry(key)
                .or_insert_with(|| Arc::new(ReadScan::from_record(&record)))
                .clone();

            let obs = scan.observation_at(qpos, base_qual);
            let bi = base_index(read_base);
            per_base[bi].push(&obs);
            depth += 1;
        }

        results.push(PositionResult {
            chrom: chrom.clone(),
            pos: one_based,
            ref_base,
            depth,
            per_base,
        });
    }

    Ok(results)
}

fn fetch_ref_base(fasta: &faidx::Reader, chrom: &str, zero_based: u32) -> Result<u8> {
    let seq = fasta
        .fetch_seq(chrom, zero_based as usize, zero_based as usize)
        .with_context(|| format!("fetching ref {}:{}", chrom, zero_based))?;
    Ok(if seq.is_empty() { b'N' } else { seq[0].to_ascii_uppercase() })
}

fn qname_hash(qname: &[u8], flag: u16) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    qname.hash(&mut h);
    flag.hash(&mut h);
    h.finish()
}

pub fn _silence_unused() {
    let _ = bam::Header::new();
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
