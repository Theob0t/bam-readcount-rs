mod bed;
mod metrics;
mod pileup;

use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::bounded;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "bam-readcount-rs",
    version,
    about = "Fast Rust reimplementation of bam-readcount, drop-in compatible with the SNV per-base output format consumed by STREGA."
)]
struct Args {
    /// Reference FASTA (must be .fai-indexed alongside).
    #[arg(short = 'f', long = "reference-fasta")]
    reference: PathBuf,

    /// Site list file (BED-like). Each line: `chr start end` where start==end
    /// for single-position queries, or a half-open range otherwise.
    #[arg(short = 'l', long = "site-list")]
    site_list: PathBuf,

    /// Minimum mapping quality.
    #[arg(short = 'q', long = "min-mapping-quality", default_value_t = 0)]
    min_mapping_quality: u8,

    /// Minimum base quality.
    #[arg(short = 'b', long = "min-base-quality", default_value_t = 0)]
    min_base_quality: u8,

    /// Max depth (clamps the pileup column).
    #[arg(short = 'd', long = "max-count", default_value_t = 10_000_000u32)]
    max_count: u32,

    /// Window size around each site. Always 1 for STREGA usage; accepted but
    /// not used (matches the `-w1` flag passed by the existing pipeline).
    #[arg(short = 'w', long = "max-warnings", default_value_t = 1)]
    _max_warnings: u32,

    /// Number of worker threads (defaults to 1; STREGA's
    /// bamreadscounts_parallel.py used 22 subprocess copies — set --threads
    /// directly here to replace that fan-out).
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    threads: usize,

    /// Output file. If omitted, writes to stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// BAM file (positional).
    bam: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .context("building rayon thread pool")?;

    let mut sites = bed::load_sites(&args.site_list)?;
    sites.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.pos.cmp(&b.pos)));
    // 200 sites per chunk, split if any gap > 50kb. Exonic BEDs cluster ~1
    // site per few hundred bp inside genes, then a multi-Mb gap to the next
    // gene; the gap rule ensures each chunk fetches a compact range.
    let chunks = bed::chunk_sites(sites, 200, 50_000);
    let n_chunks = chunks.len();

    let config = pileup::PileupConfig {
        min_mapq: args.min_mapping_quality,
        min_baseq: args.min_base_quality,
        max_depth: args.max_count,
    };

    // Each chunk emits a (chunk_index, formatted_string) so the writer thread
    // can stitch them back in BED order.
    let (tx, rx) = bounded::<(usize, String)>(n_chunks.max(1));

    let bam_path = args.bam.clone();
    let fasta_path = args.reference.clone();

    let writer_thread = std::thread::spawn(move || -> Result<()> {
        let mut out: Box<dyn Write + Send> = match args.output {
            Some(p) => Box::new(BufWriter::with_capacity(
                1 << 20,
                std::fs::File::create(&p)
                    .with_context(|| format!("creating {}", p.display()))?,
            )),
            None => Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout())),
        };

        // Reorder chunk outputs into BED order. We expect at most n_chunks
        // messages.
        let mut buffer: Vec<Option<String>> = vec![None; n_chunks];
        let mut next_to_emit: usize = 0;

        for (idx, body) in rx.iter() {
            buffer[idx] = Some(body);
            while next_to_emit < n_chunks {
                if let Some(body) = buffer[next_to_emit].take() {
                    out.write_all(body.as_bytes())?;
                    next_to_emit += 1;
                } else {
                    break;
                }
            }
        }
        out.flush()?;
        Ok(())
    });

    pool.install(|| {
        chunks
            .into_par_iter()
            .enumerate()
            .for_each(|(idx, chunk)| {
                let mut body = String::with_capacity(chunk.len() * 200);
                match pileup::process_chunk(&bam_path, &fasta_path, &chunk, config) {
                    Ok(results) => {
                        for r in results {
                            pileup::format_position(&r, &mut body);
                        }
                    }
                    Err(e) => {
                        eprintln!("[bam-readcount-rs] chunk {} failed: {:?}", idx, e);
                    }
                }
                let _ = tx.send((idx, body));
            });
    });

    drop(tx);
    writer_thread
        .join()
        .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;

    Ok(())
}
