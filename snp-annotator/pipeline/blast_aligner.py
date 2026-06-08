"""
blast_aligner.py
────────────────
Submits the Sanger sequence to NCBI BLAST (megablast against nt/human)
and parses the top hit to extract:
  - chromosome number
  - genomic coordinates (hsp sbjct_start / sbjct_end)
  - strand direction
  - list of mismatch positions with ref and query bases
"""

import time
import re
from Bio.Blast import NCBIWWW, NCBIXML


# BLAST parameters
BLAST_PROGRAM  = "blastn"
BLAST_DATABASE = "nt"          # nucleotide collection
MEGABLAST      = True
ENTREZ_QUERY   = 'Homo sapiens[Organism]'


def run_blast(sequence: str) -> dict | None:
    """
    Run NCBI BLAST and return alignment result dict.

    Returns
    -------
    {
      "chromosome":  str,           e.g. "17"
      "hsp": {
          "sbjct_start": int,
          "sbjct_end":   int,
          "query":       str,       aligned query string (with gaps)
          "sbjct":       str,       aligned subject string
          "match":       str,       match string (|, spaces, mismatches)
          "score":       float,
          "expect":      float,
          "identities":  int,
          "gaps":        int,
      },
      "strand":      int,           +1 or -1
      "mismatches":  int,           count
      "raw_mismatches": list[dict]  raw position data for dbSNP lookup
    }
    or None if no hits found.
    """
    # NCBI BLAST web API call (synchronous, waits for results)
    result_handle = NCBIWWW.qblast(
        program=BLAST_PROGRAM,
        database=BLAST_DATABASE,
        sequence=sequence,
        megablast=MEGABLAST,
        entrez_query=ENTREZ_QUERY,
        hitlist_size=5,
        alignments=5,
        descriptions=5,
    )

    blast_records = list(NCBIXML.parse(result_handle))
    if not blast_records or not blast_records[0].alignments:
        return None

    blast_record = blast_records[0]
    alignment    = blast_record.alignments[0]
    hsp          = alignment.hsps[0]

    # Parse chromosome from the hit title
    chrom = _extract_chromosome(alignment.title)

    # Determine strand
    strand = 1 if hsp.sbjct_start < hsp.sbjct_end else -1
    sbjct_start = min(hsp.sbjct_start, hsp.sbjct_end)
    sbjct_end   = max(hsp.sbjct_start, hsp.sbjct_end)

    # Extract mismatches
    raw_mismatches, mismatch_count = _extract_mismatches(
        hsp, sbjct_start, strand
    )

    return {
        "chromosome": chrom,
        "hsp": {
            "sbjct_start": sbjct_start,
            "sbjct_end":   sbjct_end,
            "query":       hsp.query,
            "sbjct":       hsp.sbjct,
            "match":       hsp.match,
            "score":       hsp.score,
            "expect":      hsp.expect,
            "identities":  hsp.identities,
            "gaps":        hsp.gaps,
        },
        "strand":          strand,
        "mismatches":      mismatch_count,
        "raw_mismatches":  raw_mismatches,
        "query_length":    blast_record.query_length,
        "accession":       alignment.accession,
    }


def _extract_chromosome(title: str) -> str:
    """
    Try to pull a chromosome number/name from a BLAST hit title.
    Examples:
      'Homo sapiens chromosome 17, GRCh38.p14 Primary Assembly'  → '17'
      'NC_000017.11 Homo sapiens chromosome 17 ...'              → '17'
    """
    # Named chromosomes
    m = re.search(r'chromosome\s+(\d+|X|Y|MT|M)', title, re.IGNORECASE)
    if m:
        return m.group(1)

    # NC_ accession  NC_000017 → chromosome 17
    m = re.search(r'NC_0000(\d{2})', title)
    if m:
        n = int(m.group(1))
        if n <= 22:
            return str(n)
        if n == 23:
            return "X"
        if n == 24:
            return "Y"
        if n == 12:
            return "MT"

    return "?"


def _extract_mismatches(hsp, sbjct_start: int, strand: int) -> tuple[list, int]:
    """
    Walk the aligned query and subject strings and collect every position
    where they differ (substitution, not a gap).

    Returns (list_of_mismatch_dicts, count)
    Each dict:
        position    : int  (genomic coordinate on subject / reference)
        query_base  : str  (base from the .ab1 sequence)
        ref_base    : str  (base from the reference genome)
        query_index : int  (0-based index in original query sequence)
    """
    mismatches = []
    query_str  = hsp.query
    sbjct_str  = hsp.sbjct
    match_str  = hsp.match

    genomic_pos  = sbjct_start
    query_offset = 0

    for i, (q, s, m) in enumerate(zip(query_str, sbjct_str, match_str)):
        if q == '-':
            # insertion in query — advance genomic, not query
            if strand == 1:
                genomic_pos += 1
            continue
        if s == '-':
            # deletion in query — advance query only
            query_offset += 1
            continue

        if m == ' ' or (q.upper() != s.upper() and m != '|'):
            # Mismatch
            mismatches.append({
                "position":    genomic_pos,
                "query_base":  q.upper(),
                "ref_base":    s.upper(),
                "query_index": hsp.query_start + query_offset - 1,
            })

        if strand == 1:
            genomic_pos += 1
        else:
            genomic_pos -= 1
        query_offset += 1

    return mismatches, len(mismatches)
