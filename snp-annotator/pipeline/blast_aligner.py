"""
blast_aligner.py
"""

import time
import re
from Bio.Blast import NCBIWWW, NCBIXML

BLAST_PROGRAM  = "blastn"
BLAST_DATABASE = "nt"
MEGABLAST      = True
ENTREZ_QUERY   = 'Homo sapiens[Organism]'


def run_blast(sequence: str) -> dict | None:
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

    chrom = _extract_chromosome(alignment.title)

    strand = 1 if hsp.sbjct_start < hsp.sbjct_end else -1
    sbjct_start = min(hsp.sbjct_start, hsp.sbjct_end)
    sbjct_end   = max(hsp.sbjct_start, hsp.sbjct_end)

    raw_mismatches, mismatch_count = _extract_mismatches(hsp, sbjct_start, strand)

    return {
        "chromosome":      chrom,
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
        "hit_title":       alignment.title,
    }


def _extract_chromosome(title: str) -> str:
    # Try explicit chromosome mention first
    m = re.search(r'chromosome\s+(\d+|X|Y|MT|M)\b', title, re.IGNORECASE)
    if m:
        return m.group(1).upper()

    # Try NC_ accession number
    m = re.search(r'NC_(\d+)\.\d+', title)
    if m:
        acc_num = int(m.group(1))
        if 1 <= acc_num <= 22:
            return str(acc_num)
        if acc_num == 23:
            return "X"
        if acc_num == 24:
            return "Y"
        if acc_num == 12920:
            return "MT"

    return "?"


def _extract_mismatches(hsp, sbjct_start: int, strand: int) -> tuple[list, int]:
    mismatches = []
    query_str  = hsp.query
    sbjct_str  = hsp.sbjct
    match_str  = hsp.match

    genomic_pos  = sbjct_start
    query_offset = 0

    for i, (q, s, m) in enumerate(zip(query_str, sbjct_str, match_str)):
        if q == '-':
            if strand == 1:
                genomic_pos += 1
            continue
        if s == '-':
            query_offset += 1
            continue

        if m == ' ' or (q.upper() != s.upper() and m != '|'):
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
