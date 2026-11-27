//! Streaming MAF v1 stanza parser and two-row MAF block utilities.

/// A parsed `s`-line row from a MAF v1 alignment block.
#[derive(Debug, Clone)]
pub struct MafRow {
    /// Sequence/contig identifier.
    pub src: String,
    /// 0-based start offset on `strand`.
    pub start: usize,
    /// Strand (`'+'` or `'-'`).
    pub strand: char,
    /// Total ungapped length of the source sequence.
    pub src_size: usize,
    /// Gapped alignment row string.
    pub seq: String,
}

impl MafRow {
    /// Convert a 0-based ungapped index `nongap_k` along this row to a forward-strand coordinate.
    #[must_use]
    pub fn fwd_pos(&self, nongap_k: usize) -> usize {
        let strand =
            crate::coord::Strand::from_char(self.strand).unwrap_or(crate::coord::Strand::Plus);
        crate::coord::Interval::ucsc_offset_to_fwd_pos(self.start, nongap_k, strand, self.src_size)
    }
}

/// Parse a single trimmed MAF `s` line into `(MafRow, ungapped_size)`.
#[must_use]
pub fn parse_maf_s_line(trimmed: &str) -> Option<(MafRow, usize)> {
    let rest = trimmed
        .strip_prefix('s')
        .filter(|r| r.starts_with(char::is_whitespace))?;
    let mut parts = rest.split_whitespace();
    let (
        Some(src),
        Some(start_str),
        Some(size_str),
        Some(strand_str),
        Some(src_size_str),
        Some(seq),
    ) = (
        parts.next(),
        parts.next(),
        parts.next(),
        parts.next(),
        parts.next(),
        parts.next(),
    )
    else {
        return None;
    };
    let (Ok(start), Ok(size), Ok(src_size)) = (
        start_str.parse::<usize>(),
        size_str.parse::<usize>(),
        src_size_str.parse::<usize>(),
    ) else {
        return None;
    };
    let strand = match strand_str {
        "+" => '+',
        "-" => '-',
        _ => return None,
    };
    if start.checked_add(size).is_none_or(|end| end > src_size) {
        return None;
    }
    if !seq.is_ascii() || seq.bytes().filter(|&b| b != b'-').count() != size {
        return None;
    }
    Some((
        MafRow {
            src: src.to_owned(),
            start,
            strand,
            src_size,
            seq: seq.to_owned(),
        },
        size,
    ))
}

/// Maximum single-line byte length accepted by [`for_each_maf_stanza_reader`] (256 MiB).
pub(crate) const MAX_MAF_LINE_BYTES: usize = 256 * 1024 * 1024;

fn feed_maf_line(
    line: &str,
    rows: &mut Vec<(MafRow, usize)>,
    stanza_valid: &mut bool,
    on_stanza: &mut impl FnMut(&mut Vec<(MafRow, usize)>),
) {
    let trimmed = line.trim();
    if trimmed.starts_with('#') {
        return;
    }
    if trimmed.is_empty()
        || trimmed == "a"
        || trimmed
            .strip_prefix('a')
            .is_some_and(|r| r.starts_with(char::is_whitespace))
    {
        if *stanza_valid && !rows.is_empty() {
            on_stanza(rows);
        }
        rows.clear();
        *stanza_valid = true;
        return;
    }
    if trimmed
        .strip_prefix('s')
        .is_some_and(|r| r.starts_with(char::is_whitespace))
    {
        if let Some(parsed) = parse_maf_s_line(trimmed) {
            if *stanza_valid {
                rows.push(parsed);
            }
        } else {
            *stanza_valid = false;
        }
    }
}

/// Iterate over `a`-stanzas in a buffered MAF v1 reader, invoking `on_stanza` with parsed `(MafRow, ungapped_size)` rows.
///
/// # Errors
/// Returns an error if reading fails or a single line exceeds 256 MiB.
pub fn for_each_maf_stanza_reader<R: std::io::BufRead>(
    mut reader: R,
    mut on_stanza: impl FnMut(&mut Vec<(MafRow, usize)>),
) -> anyhow::Result<()> {
    let mut rows: Vec<(MafRow, usize)> = Vec::new();
    let mut stanza_valid = true;
    let mut line = String::new();
    let mut line_no = 0usize;
    loop {
        line.clear();
        let n = std::io::BufRead::read_line(
            &mut std::io::Read::take(&mut reader, MAX_MAF_LINE_BYTES as u64 + 1),
            &mut line,
        )?;
        if n == 0 {
            break;
        }
        line_no += 1;
        if line.len() > MAX_MAF_LINE_BYTES {
            anyhow::bail!("MAF line exceeds 256 MiB safety limit");
        }
        let trimmed = line.trim();
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            let valid_prefix = trimmed == "a"
                || trimmed
                    .strip_prefix(['a', 's', 'i', 'e', 'q'])
                    .is_some_and(|r| r.starts_with(char::is_whitespace));
            if !valid_prefix {
                anyhow::bail!("malformed MAF line {line_no}: unrecognized record prefix");
            }
        }
        feed_maf_line(&line, &mut rows, &mut stanza_valid, &mut on_stanza);
    }
    if stanza_valid && !rows.is_empty() {
        on_stanza(&mut rows);
    }
    Ok(())
}

/// Iterate over `a`-stanzas in a MAF v1 text stream, invoking `on_stanza` with the parsed `(MafRow, ungapped_size)` rows.
pub(crate) fn for_each_maf_stanza(
    maf_text: &str,
    mut on_stanza: impl FnMut(&mut Vec<(MafRow, usize)>),
) {
    let mut rows: Vec<(MafRow, usize)> = Vec::new();
    let mut stanza_valid = true;
    for line in maf_text.lines() {
        feed_maf_line(line, &mut rows, &mut stanza_valid, &mut on_stanza);
    }
    if stanza_valid && !rows.is_empty() {
        on_stanza(&mut rows);
    }
}

/// Parse two-row pairwise alignment blocks from a MAF v1 text stream.
#[must_use]
pub fn parse_two_row_maf_blocks(maf_text: &str) -> Vec<(MafRow, MafRow)> {
    let mut blocks = Vec::new();
    for_each_maf_stanza(maf_text, |rows| {
        if rows.len() >= 2 {
            let mut drain = rows.drain(..);
            if let (Some((r1, _)), Some((r2, _))) = (drain.next(), drain.next()) {
                blocks.push((r1, r2));
            }
        }
    });
    blocks
}
