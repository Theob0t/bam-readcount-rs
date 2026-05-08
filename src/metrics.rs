use rust_htslib::bam::record::Record;
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
        if let Some(v) = obs.se_mapq {
            self.sum_se_mq += v as u64;
        }
        // Else: read is properly-paired but missing SM tag — upstream skips
        // the SE-MAPQ contribution but still counts the read in `count`,
        // which drags the average DOWN. Match that exactly.
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
    /// Single-ended mapping quality contribution. None means "skip" — upstream
    /// semantics: a properly-paired read with no SM tag does NOT contribute
    /// to the running sum (but still counts toward read_count, dragging the
    /// average down).
    pub se_mapq: Option<u8>,
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
/// Computed once per read, then reused for every wanted aligned position the
/// read covers. The previous implementation recomputed this per-pileup-hit,
/// which on deep samples ran 5-10x more often than necessary.
#[derive(Debug, Clone)]
pub struct ReadScan {
    pub mapq: u8,
    pub se_mapq: Option<u8>,
    pub is_reverse: bool,
    pub l_qseq: u32,
    pub clipped_length: u32, // l_qseq - left_soft - right_soft
    pub left_clip: u32,      // = left_soft
    pub q2_pos: i32,
    pub three_prime_index: i32,
    pub nm: u32,
    /// Per-read sum of mismatch qualities. Adjacent mismatches collapse to the
    /// MAX base qual within a run (upstream behavior); reference-N (MD says
    /// '0' or skipped) and read-`=` (rare) do not contribute.
    pub sum_mismatch_quals: u32,
}

// CIGAR op codes per BAM spec.
const OP_MATCH: u8 = 0;
const OP_INS: u8 = 1;
const OP_DEL: u8 = 2;
const OP_REFSKIP: u8 = 3;
const OP_SOFT_CLIP: u8 = 4;
const OP_HARD_CLIP: u8 = 5;
#[allow(dead_code)]
const OP_PAD: u8 = 6;
const OP_EQUAL: u8 = 7;
const OP_DIFF: u8 = 8;

#[inline(always)]
pub fn cigar_decode(packed: u32) -> (u8, u32) {
    ((packed & 0xf) as u8, packed >> 4)
}

impl ReadScan {
    pub fn from_record(record: &Record) -> Self {
        let mapq = record.mapq();
        let is_reverse = record.is_reverse();
        let l_qseq = record.seq_len() as u32;

        // Soft-clip layout from raw_cigar (no Vec<Cigar> allocation).
        let raw = record.raw_cigar();
        let n = raw.len();
        let mut left_soft: u32 = 0;
        let mut right_soft: u32 = 0;
        if n >= 1 {
            let (op0, len0) = cigar_decode(raw[0]);
            if op0 == OP_SOFT_CLIP {
                left_soft += len0;
            }
            if n >= 2 && op0 == OP_HARD_CLIP {
                let (op1, len1) = cigar_decode(raw[1]);
                if op1 == OP_SOFT_CLIP {
                    left_soft += len1;
                }
            }
            let (opn, lenn) = cigar_decode(raw[n - 1]);
            if opn == OP_SOFT_CLIP {
                right_soft += lenn;
            }
            if n >= 2 && opn == OP_HARD_CLIP {
                let (opnm1, lennm1) = cigar_decode(raw[n - 2]);
                if opnm1 == OP_SOFT_CLIP {
                    right_soft += lennm1;
                }
            }
        }
        let clipped_length = l_qseq.saturating_sub(left_soft + right_soft);
        let left_clip = left_soft as i32;
        let right_clip = (l_qseq.saturating_sub(right_soft)) as i32;

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
        let q2_pos = find_q2_pos(qual, is_reverse);
        let three_prime_index =
            compute_three_prime_index(l_qseq, is_reverse, left_clip, right_clip, q2_pos);

        // MD walk — synchronized with raw_cigar in a single pass; no
        // aligned_to_read materialization.
        let sum_mismatch_quals = match record.aux(b"MD") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => {
                scan_mismatch_qualities_streaming(s.as_bytes(), qual, raw)
            }
            _ => 0,
        };

        let is_proper_pair = (record.flags() & 0x2) != 0;
        let se_mapq = if is_proper_pair {
            extract_sm_tag(record) // None when SM missing → upstream skips
        } else {
            Some(mapq)
        };

        Self {
            mapq,
            se_mapq,
            is_reverse,
            l_qseq,
            clipped_length,
            left_clip: left_soft,
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

/// Walk MD synchronously with raw cigar to compute per-read sum of mismatch
/// base qualities. Adjacent mismatches (consecutive read positions) collapse
/// to MAX qual within a run (upstream bamreadcount.cpp:138-199 behavior).
///
/// Streaming: no `aligned_to_read: Vec<usize>` allocation, no MD pre-parse. We
/// keep a tiny state machine over the MD bytes (number-of-pending-matches +
/// byte cursor) and step it by 1 each time we encounter an M/=/X aligned
/// position in CIGAR.
fn scan_mismatch_qualities_streaming(md: &[u8], qual: &[u8], raw_cigar: &[u32]) -> u32 {
    let mut sum: u32 = 0;
    let mut last_mm_pos: i32 = -1;
    let mut last_mm_qual: u32 = 0;

    let mut md_pos: usize = 0;
    let mut pending_match: u32 = 0;
    let mut read_idx: u32 = 0;

    for &packed in raw_cigar {
        let (op, len) = cigar_decode(packed);
        match op {
            OP_MATCH | OP_EQUAL | OP_DIFF => {
                for _ in 0..len {
                    let is_mismatch = md_advance_one(md, &mut md_pos, &mut pending_match);
                    if is_mismatch {
                        if (read_idx as usize) < qual.len() {
                            let q = qual[read_idx as usize] as u32;
                            if last_mm_pos != -1 && (last_mm_pos as u32) + 1 == read_idx {
                                if q > last_mm_qual {
                                    last_mm_qual = q;
                                }
                            } else {
                                if last_mm_pos != -1 {
                                    sum += last_mm_qual;
                                }
                                last_mm_qual = q;
                            }
                            last_mm_pos = read_idx as i32;
                        }
                    }
                    read_idx += 1;
                }
            }
            OP_INS | OP_SOFT_CLIP => {
                read_idx += len;
            }
            OP_DEL => {
                md_consume_deletion(md, &mut md_pos, len);
            }
            // RefSkip/HardClip/Pad: no MD or read advance for our purposes.
            _ => {}
        }
    }

    if last_mm_pos != -1 {
        sum += last_mm_qual;
    }
    sum
}

#[inline]
fn md_advance_one(md: &[u8], pos: &mut usize, pending_match: &mut u32) -> bool {
    if *pending_match > 0 {
        *pending_match -= 1;
        return false;
    }
    while *pos < md.len() {
        let c = md[*pos];
        if c.is_ascii_digit() {
            let mut n: u32 = 0;
            while *pos < md.len() && md[*pos].is_ascii_digit() {
                n = n.saturating_mul(10).saturating_add((md[*pos] - b'0') as u32);
                *pos += 1;
            }
            if n == 0 {
                continue; // zero-length match, retry
            }
            *pending_match = n - 1;
            return false;
        } else if c == b'^' {
            // Defensive: shouldn't appear here (caller should consume_deletion
            // before walking past a CIGAR D op). Skip silently.
            *pos += 1;
            while *pos < md.len() && md[*pos].is_ascii_alphabetic() {
                *pos += 1;
            }
            continue;
        } else if c.is_ascii_alphabetic() {
            *pos += 1;
            return true;
        } else {
            *pos += 1;
        }
    }
    false // MD exhausted — treat as match
}

#[inline]
fn md_consume_deletion(md: &[u8], pos: &mut usize, len: u32) {
    while *pos < md.len() && md[*pos] == b'0' {
        *pos += 1;
    }
    if *pos < md.len() && md[*pos] == b'^' {
        *pos += 1;
        let mut consumed: u32 = 0;
        while consumed < len && *pos < md.len() && md[*pos].is_ascii_alphabetic() {
            *pos += 1;
            consumed += 1;
        }
    }
}

/// Find q2_pos exactly as upstream's fetch_func does. For reverse reads
/// the scan starts at `k = 0` (NOT `left_clip`) — upstream walks through
/// the soft-clipped region too, so a non-Q2 base at qpos 0 yields
/// `q2_pos = -1`.
fn find_q2_pos(qual: &[u8], is_reverse: bool) -> i32 {
    let l = qual.len() as i32;
    if l == 0 {
        return -1;
    }
    let mut q2_pos: i32 = -1;
    let (mut k, increment): (i32, i32) = if is_reverse { (0, 1) } else { (l - 1, -1) };
    while q2_pos < 0 && k >= 0 && k < l {
        if qual[k as usize] != 2 {
            q2_pos = k - 1;
            break;
        }
        k += increment;
    }
    q2_pos
}

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

fn extract_sm_tag(record: &Record) -> Option<u8> {
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
