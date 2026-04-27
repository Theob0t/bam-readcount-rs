use rust_htslib::bam::record::{Cigar, Record};
use std::fmt::Write;

/// One per-base record at one position. The 13 metrics after `count` follow the
/// upstream bam-readcount v1.0.1 ordering, matching `posLevel.read_bamReadCountsFile`
/// (`STREGA/STREGA/posLevel.py:207`).
///
/// Formulas reproduce upstream `BasicStat::process_read` and `output_stat` in
/// `src/lib/bamrc/BasicStat.cpp`, with the per-read Zm-tag values computed as
/// in `src/exe/bam-readcount/bamreadcount.cpp::fetch_func()`.
#[derive(Debug, Clone, Default)]
pub struct PerBaseRecord {
    pub count: u32,
    sum_mq: u64,
    sum_bq: u64,
    sum_se_mq: u64,
    pub num_plus: u32,
    pub num_minus: u32,
    sum_pos_frac: f64,
    sum_mm_frac: f64,
    sum_sum_mm_q: f64,
    pub num_q2_reads: u32,
    sum_dist_q2: f64, // averaged only over reads with q2_val > -1
    sum_clipped_len: u64,
    sum_dist_3p: f64,
}

impl PerBaseRecord {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn push(&mut self, obs: &ReadObservation) {
        self.count += 1;
        self.sum_mq += obs.mapq as u64;
        self.sum_bq += obs.base_qual as u64;
        self.sum_se_mq += obs.se_mapq as u64;
        if obs.is_reverse {
            self.num_minus += 1;
        } else {
            self.num_plus += 1;
        }
        self.sum_pos_frac += obs.pos_as_fraction;
        self.sum_mm_frac += obs.mismatch_fraction;
        self.sum_sum_mm_q += obs.sum_mismatch_quals as f64;
        if let Some(d) = obs.q2_distance {
            self.num_q2_reads += 1;
            self.sum_dist_q2 += d;
        }
        self.sum_clipped_len += obs.clipped_len as u64;
        self.sum_dist_3p += obs.distance_to_3p_end;
    }

    pub fn write_with_label(&self, base: &str, out: &mut String) {
        if self.count == 0 {
            out.push_str(base);
            out.push_str(":0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00");
            return;
        }
        let n = self.count as f64;
        let avg_mq = self.sum_mq as f64 / n;
        let avg_bq = self.sum_bq as f64 / n;
        let avg_se_mq = self.sum_se_mq as f64 / n;
        let avg_pos_frac = self.sum_pos_frac / n;
        let avg_mm_frac = self.sum_mm_frac / n;
        let avg_sum_mm_q = self.sum_sum_mm_q / n;
        let avg_dist_q2 = if self.num_q2_reads > 0 {
            self.sum_dist_q2 / self.num_q2_reads as f64
        } else {
            0.0
        };
        let avg_clipped = self.sum_clipped_len as f64 / n;
        let avg_dist_3p = self.sum_dist_3p / n;

        write!(
            out,
            "{}:{}:{:.2}:{:.2}:{:.2}:{}:{}:{:.2}:{:.2}:{:.2}:{}:{:.2}:{:.2}:{:.2}",
            base,
            self.count,
            avg_mq,
            avg_bq,
            avg_se_mq,
            self.num_plus,
            self.num_minus,
            avg_pos_frac,
            avg_mm_frac,
            avg_sum_mm_q,
            self.num_q2_reads,
            avg_dist_q2,
            avg_clipped,
            avg_dist_3p
        )
        .expect("writing to String never fails");
    }
}

/// Per-position observation for one read, ready to push into a PerBaseRecord.
#[derive(Debug, Clone)]
pub struct ReadObservation {
    pub mapq: u8,
    pub se_mapq: u8,
    pub base_qual: u8,
    pub is_reverse: bool,
    pub pos_as_fraction: f64,
    pub mismatch_fraction: f64,
    pub sum_mismatch_quals: u32,
    pub q2_distance: Option<f64>,
    pub clipped_len: u32, // = `clipped_length` upstream (aligned, soft-clip-removed)
    pub distance_to_3p_end: f64,
}

/// Per-read scratch values (the upstream "Zm" tag fields, computed inline).
/// Computed once per read regardless of how many of the read's aligned bases
/// land in the BED.
#[derive(Debug, Clone)]
pub struct ReadScan {
    pub mapq: u8,
    pub se_mapq: u8,
    pub is_reverse: bool,
    pub l_qseq: u32,
    pub clipped_length: u32, // l_qseq - left_soft - right_soft
    pub left_clip: u32,      // = left_soft
    #[allow(dead_code)]
    pub right_clip: u32,     // upstream coord = l_qseq - right_soft (kept for symmetry/debug)
    /// q2_pos in qpos coordinates (signed; -1 means none). For forward reads
    /// `q2_pos = l_qseq - 2` whenever the very last base is not Q2 (this is
    /// the upstream behavior, see fetch_func() in bamreadcount.cpp).
    pub q2_pos: i32,
    pub three_prime_index: i32,
    pub nm: u32,
    /// Per-read sum of mismatch qualities. Adjacent mismatches collapse to the
    /// MAX base qual within a run (upstream behavior); reference-N (MD says
    /// '0' or skipped) and read-`=` (rare) do not contribute.
    pub sum_mismatch_quals: u32,
}

impl ReadScan {
    pub fn from_record(record: &Record) -> Self {
        let mapq = record.mapq();
        let is_reverse = record.is_reverse();
        let l_qseq = record.seq_len() as u32;

        // Walk CIGAR for left_soft / right_soft. Hard clips are NOT subtracted
        // (and don't appear in l_qseq).
        let cigar = record.cigar();
        let cigar_view: Vec<Cigar> = cigar.iter().copied().collect();
        let mut left_soft: u32 = 0;
        let mut right_soft: u32 = 0;
        let n = cigar_view.len();
        for (i, op) in cigar_view.iter().enumerate() {
            if let Cigar::SoftClip(l) = op {
                if i == 0 || (i == 1 && matches!(cigar_view[0], Cigar::HardClip(_))) {
                    left_soft += *l;
                } else if i == n - 1
                    || (i == n - 2 && matches!(cigar_view[n - 1], Cigar::HardClip(_)))
                {
                    right_soft += *l;
                }
            }
        }
        let clipped_length = l_qseq.saturating_sub(left_soft + right_soft);
        let left_clip = left_soft;
        let right_clip = l_qseq.saturating_sub(right_soft);

        // NM tag — used only for avg_num_mismatches_as_fraction.
        let nm = match record.aux(b"NM") {
            Ok(rust_htslib::bam::record::Aux::I32(n)) => n as u32,
            Ok(rust_htslib::bam::record::Aux::U32(n)) => n,
            Ok(rust_htslib::bam::record::Aux::I16(n)) => n.max(0) as u32,
            Ok(rust_htslib::bam::record::Aux::U16(n)) => n as u32,
            Ok(rust_htslib::bam::record::Aux::I8(n)) => n.max(0) as u32,
            Ok(rust_htslib::bam::record::Aux::U8(n)) => n as u32,
            _ => 0,
        };

        let qual = record.qual();
        let q2_pos = find_q2_pos(qual, is_reverse, left_clip);
        let three_prime_index = compute_three_prime_index(l_qseq, is_reverse, left_clip as i32, right_clip as i32, q2_pos);

        let sum_mismatch_quals = scan_mismatch_qualities_via_md(record, qual, &cigar_view);

        Self {
            mapq,
            se_mapq: extract_se_mapq(record).unwrap_or(mapq),
            is_reverse,
            l_qseq,
            clipped_length,
            left_clip,
            right_clip,
            q2_pos,
            three_prime_index,
            nm,
            sum_mismatch_quals,
        }
    }

    pub fn observation_at(&self, qpos: u32, base_qual: u8) -> ReadObservation {
        let l_qseq = self.l_qseq.max(1) as f64;
        // pf: 1 - |((qpos - left_clip) - center) / center| with center = clipped/2
        let pos_as_fraction = if self.clipped_length > 0 {
            let center = self.clipped_length as f64 / 2.0;
            let aligned_pos = qpos as i64 - self.left_clip as i64;
            1.0 - ((aligned_pos as f64 - center).abs() / center)
        } else {
            0.0
        };

        // Mismatch fraction: NM / clipped_length (upstream uses the NM aux tag)
        let mismatch_fraction = if self.clipped_length > 0 {
            self.nm as f64 / self.clipped_length as f64
        } else {
            0.0
        };

        // Q2 distance: |qpos - q2_pos| / l_qseq, only if q2_pos > -1.
        let q2_distance = if self.q2_pos > -1 {
            Some(((qpos as i64 - self.q2_pos as i64).abs() as f64) / l_qseq)
        } else {
            None
        };

        // 3p distance: |qpos - three_prime_index| / l_qseq.
        let distance_to_3p_end =
            ((qpos as i64 - self.three_prime_index as i64).abs() as f64) / l_qseq;

        ReadObservation {
            mapq: self.mapq,
            se_mapq: self.se_mapq,
            base_qual,
            is_reverse: self.is_reverse,
            pos_as_fraction,
            mismatch_fraction,
            sum_mismatch_quals: self.sum_mismatch_quals,
            q2_distance,
            clipped_len: self.clipped_length,
            distance_to_3p_end,
        }
    }
}

/// Walk MD tag with CIGAR + qualities. Returns the per-read sum of mismatch
/// base qualities, where adjacent mismatches collapse to MAX base qual within
/// a run (upstream bamreadcount.cpp:138-199 behavior).
fn scan_mismatch_qualities_via_md(record: &Record, qual: &[u8], cigar: &[Cigar]) -> u32 {
    let md = match record.aux(b"MD") {
        Ok(rust_htslib::bam::record::Aux::String(s)) => s.to_string(),
        _ => return 0,
    };

    // Build mapping aligned_idx -> read_idx (BAM seq index). MATCH/EQUAL/DIFF
    // ops consume both. INS/SOFTCLIP consume read only. DEL/REF_SKIP consume ref only.
    let seq_len = record.seq_len();
    let mut aligned_to_read: Vec<usize> = Vec::with_capacity(seq_len);
    let mut read_idx: usize = 0;
    for op in cigar {
        match op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                for _ in 0..*l {
                    aligned_to_read.push(read_idx);
                    read_idx += 1;
                }
            }
            Cigar::Ins(l) | Cigar::SoftClip(l) => {
                read_idx += *l as usize;
            }
            _ => {}
        }
    }

    let mut sum: u32 = 0;
    let mut last_mm_pos: i32 = -1;
    let mut last_mm_qual: u32 = 0;

    let bytes = md.as_bytes();
    let mut aligned_cursor: usize = 0;
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        if c.is_ascii_digit() {
            let mut j = i;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            let n: usize = std::str::from_utf8(&bytes[i..j]).unwrap_or("0").parse().unwrap_or(0);
            aligned_cursor += n;
            i = j;
        } else if c == b'^' {
            i += 1;
            while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                i += 1;
            }
        } else if c.is_ascii_alphabetic() {
            // Mismatch at current aligned position. Look up the read position.
            if aligned_cursor < aligned_to_read.len() {
                let read_pos = aligned_to_read[aligned_cursor];
                if read_pos < qual.len() {
                    let q = qual[read_pos] as u32;
                    if last_mm_pos != -1 {
                        if last_mm_pos as usize + 1 != read_pos {
                            // previous run ends; flush
                            sum += last_mm_qual;
                            last_mm_qual = q;
                            last_mm_pos = read_pos as i32;
                        } else {
                            // adjacent: keep the maximum qual in this run
                            if last_mm_qual < q {
                                last_mm_qual = q;
                            }
                            last_mm_pos = read_pos as i32;
                        }
                    } else {
                        last_mm_pos = read_pos as i32;
                        last_mm_qual = q;
                    }
                }
            }
            aligned_cursor += 1;
            i += 1;
        } else {
            i += 1;
        }
    }
    if last_mm_pos != -1 {
        sum += last_mm_qual;
    }
    sum
}

/// Find q2_pos exactly as upstream's fetch_func does.
fn find_q2_pos(qual: &[u8], is_reverse: bool, left_clip: u32) -> i32 {
    let l = qual.len() as i32;
    if l == 0 {
        return -1;
    }
    let mut q2_pos: i32 = -1;
    let (mut k, increment): (i32, i32) = if is_reverse {
        (left_clip as i32, 1)
    } else {
        (l - 1, -1)
    };
    while q2_pos < 0 && k >= 0 && k < l {
        if qual[k as usize] != 2 {
            q2_pos = k - 1;
            break;
        }
        k += increment;
    }
    q2_pos
}

/// `three_prime_index` matches upstream bamreadcount.cpp: clamped to right_clip
/// (forward) or left_clip (reverse), then pulled inward by q2_pos.
fn compute_three_prime_index(
    l_qseq: u32,
    is_reverse: bool,
    left_clip: i32,
    right_clip: i32,
    q2_pos: i32,
) -> i32 {
    let mut tpi: i32;
    if is_reverse {
        tpi = 0;
        if tpi < left_clip {
            tpi = left_clip;
        }
        if q2_pos > tpi {
            tpi = q2_pos;
        }
    } else {
        tpi = (l_qseq as i32) - 1;
        if tpi > right_clip {
            tpi = right_clip;
        }
        if tpi > q2_pos && q2_pos != -1 {
            tpi = q2_pos;
        }
    }
    tpi
}

fn extract_se_mapq(record: &Record) -> Option<u8> {
    // Upstream uses the "SM" (single-ended mapping quality) tag if present,
    // else falls back to MAPQ. Note this is "SM" not "SA" — a different tag.
    match record.aux(b"SM") {
        Ok(rust_htslib::bam::record::Aux::I32(n)) => Some(n.clamp(0, 255) as u8),
        Ok(rust_htslib::bam::record::Aux::U32(n)) => Some(n.min(255) as u8),
        Ok(rust_htslib::bam::record::Aux::I16(n)) => Some(n.clamp(0, 255) as u8),
        Ok(rust_htslib::bam::record::Aux::U16(n)) => Some(n.min(255) as u8),
        Ok(rust_htslib::bam::record::Aux::I8(n)) => Some(n.max(0) as u8),
        Ok(rust_htslib::bam::record::Aux::U8(n)) => Some(n),
        _ => None,
    }
}

pub fn base_index(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

pub const BASE_LABELS: [&str; 5] = ["A", "C", "G", "T", "N"];
