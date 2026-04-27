use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Site {
    pub chrom: String,
    pub pos: u32,
}

pub fn load_sites<P: AsRef<Path>>(path: P) -> Result<Vec<Site>> {
    let file = File::open(path.as_ref())
        .with_context(|| format!("opening BED {}", path.as_ref().display()))?;
    let reader = BufReader::new(file);

    let mut sites = Vec::new();
    for (i, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let mut fields = line.split('\t');
        let chrom = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("BED line {} missing chrom", i + 1))?
            .to_string();
        let start: u32 = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("BED line {} missing start", i + 1))?
            .parse()
            .with_context(|| format!("parsing start on BED line {}", i + 1))?;
        let end: u32 = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("BED line {} missing end", i + 1))?
            .parse()
            .with_context(|| format!("parsing end on BED line {}", i + 1))?;

        // bam-readcount semantics: site list lines specify single positions.
        // STREGA emits 1-based "chr pos pos" (start == end). We emit one site
        // per BED line, at the "start" (which for STREGA is the 1-based pos).
        // For genuine BED ranges (start..end half-open), expand to (end-start)
        // sites using start+1..=end (1-based) when start < end.
        if start == end {
            sites.push(Site { chrom, pos: start });
        } else {
            for p in (start + 1)..=end {
                sites.push(Site {
                    chrom: chrom.clone(),
                    pos: p,
                });
            }
        }
    }
    Ok(sites)
}

/// Group sites into chunks suitable for one BAM fetch + pileup pass each.
///
/// A new chunk starts when any of:
///   - chrom changes
///   - chunk already has `target_size` positions
///   - gap from previous site exceeds `max_gap` bp (so each chunk has a
///     compact genomic span; pileup walks every column in the fetch range,
///     so a 100-position chunk spread over 100 Mbp would be catastrophically
///     slow)
///
/// Sites must already be sorted by (chrom, pos).
pub fn chunk_sites(sites: Vec<Site>, target_size: usize, max_gap: u32) -> Vec<Vec<Site>> {
    if sites.is_empty() {
        return Vec::new();
    }
    let mut chunks: Vec<Vec<Site>> = Vec::new();
    let mut current: Vec<Site> = Vec::with_capacity(target_size);
    let mut current_chrom = sites[0].chrom.clone();
    let mut last_pos: u32 = 0;

    for site in sites {
        let split = site.chrom != current_chrom
            || current.len() >= target_size
            || (!current.is_empty() && site.pos.saturating_sub(last_pos) > max_gap);
        if split {
            if !current.is_empty() {
                chunks.push(std::mem::take(&mut current));
            }
            current_chrom = site.chrom.clone();
        }
        last_pos = site.pos;
        current.push(site);
    }
    if !current.is_empty() {
        chunks.push(current);
    }
    chunks
}
