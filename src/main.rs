mod bed;
mod metrics;
mod pileup;

use anyhow::{Context, Result};
use clap::Parser;
use crossbeam_channel::bounded;
use rust_htslib::bam::{IndexedReader, Read as BamRead};
use rust_htslib::faidx;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::bed::Site;
use crate::pileup::ChunkScratch;

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

    /// Total CPU budget. Workers × bgzf-threads should equal this.
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    threads: usize,

    /// htslib BGZF decompression threads per worker. Multi-threaded BGZF is
    /// often the single biggest win on deep BAMs (BAM I/O is the bottleneck
    /// once ReadScan stops being recomputed). With `--threads N` and
    /// `--bgzf-threads K`, we spawn `max(1, N / K)` worker threads, each
    /// owning one IndexedReader configured with K background BGZF threads.
    #[arg(long = "bgzf-threads", default_value_t = 2)]
    bgzf_threads: usize,

    /// Output file. If omitted, writes to stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// BAM file (positional).
    bam: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();

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

    // Worker pool: each worker keeps one IndexedReader + one faidx + one
    // ChunkScratch open for the entire run, so we don't pay BAM/index reopen
    // cost per chunk. With persistent BGZF threads via set_threads(K), this
    // is the configuration that wins on deep samples.
    let bgzf_threads = args.bgzf_threads.max(1);
    let n_workers = (args.threads / bgzf_threads).max(1);

    let (work_tx, work_rx) = bounded::<(usize, Vec<Site>)>(n_chunks.max(1));
    let (result_tx, result_rx) = bounded::<(usize, String)>(n_chunks.max(1));

    let bam_path = args.bam.clone();
    let fasta_path = args.reference.clone();

    let mut workers = Vec::with_capacity(n_workers);
    for _ in 0..n_workers {
        let work_rx = work_rx.clone();
        let result_tx = result_tx.clone();
        let bam_path = bam_path.clone();
        let fasta_path = fasta_path.clone();
        let bgzf = bgzf_threads;
        let h = std::thread::spawn(move || -> Result<()> {
            let mut reader = IndexedReader::from_path(&bam_path)
                .with_context(|| format!("opening BAM {}", bam_path.display()))?;
            // Single-threaded set_threads(0) is a no-op; >0 enables htslib's
            // BGZF thread pool. Only call when >1 to keep behavior identical
            // to the old single-threaded path in 1-cpu runs.
            if bgzf > 1 {
                reader
                    .set_threads(bgzf)
                    .context("enabling BGZF threads on IndexedReader")?;
            }
            let fasta = faidx::Reader::from_path(&fasta_path)
                .with_context(|| format!("opening FASTA {}", fasta_path.display()))?;
            let mut scratch = ChunkScratch::new();

            for (idx, chunk) in work_rx.iter() {
                let mut body = String::with_capacity(chunk.len() * 200);
                match pileup::process_chunk(
                    &mut reader,
                    &fasta,
                    &chunk,
                    config,
                    &mut scratch,
                ) {
                    Ok(results) => {
                        for r in results {
                            pileup::format_position(&r, &mut body);
                        }
                    }
                    Err(e) => {
                        eprintln!("[bam-readcount-rs] chunk {} failed: {:?}", idx, e);
                    }
                }
                if result_tx.send((idx, body)).is_err() {
                    break; // writer thread closed
                }
            }
            Ok(())
        });
        workers.push(h);
    }

    // Drop our handles so the channels close once all workers finish.
    drop(work_rx);
    drop(result_tx);

    // Push work.
    for (idx, chunk) in chunks.into_iter().enumerate() {
        if work_tx.send((idx, chunk)).is_err() {
            break;
        }
    }
    drop(work_tx);

    // Writer thread: reorders chunk outputs back into BED order.
    let output = args.output;
    let writer_thread = std::thread::spawn(move || -> Result<()> {
        let mut out: Box<dyn Write + Send> = match output {
            Some(p) => Box::new(BufWriter::with_capacity(
                1 << 20,
                std::fs::File::create(&p)
                    .with_context(|| format!("creating {}", p.display()))?,
            )),
            None => Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout())),
        };

        let mut buffer: Vec<Option<String>> = vec![None; n_chunks];
        let mut next_to_emit: usize = 0;

        for (idx, body) in result_rx.iter() {
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

    for h in workers {
        h.join()
            .map_err(|_| anyhow::anyhow!("worker thread panicked"))??;
    }
    writer_thread
        .join()
        .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;

    Ok(())
}
