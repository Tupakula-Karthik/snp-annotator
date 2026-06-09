"""
blast_aligner.py — uses NCBI BLAST REST API instead of BioPython qblast
"""

import time
import re
import json
import urllib.request
import urllib.parse

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"


def run_blast(sequence: str) -> dict | None:
    sequence = sequence.strip("N")
    if len(sequence) < 100:
        raise ValueError("Sequence too short (need 100+ bp)")

    # Step 1 — Submit job
    rid = _submit_blast(sequence)
    if not rid:
        raise ValueError("Failed to submit BLAST job")

    # Step 2 — Poll until done
    results_xml = _poll_blast(rid)
    if not results_xml:
        raise ValueError("BLAST timed out or returned no results")

    # Step 3 — Parse XML
    return _parse_blast_xml(results_xml)


def _submit_blast(sequence: str) -> str | None:
    params = urllib.parse.urlencode({
        "CMD":        "Put",
        "PROGRAM":    "blastn",
        "DATABASE":   "nt",
        "QUERY":      sequence,
        "MEGABLAST":  "on",
        "ENTREZ_QUERY": "Homo sapiens[Organism]",
        "HITLIST_SIZE": "5",
        "FORMAT_TYPE":  "XML",
    }).encode()

    req = urllib.request.Request(BLAST_URL, data=params)
    with urllib.request.urlopen(req, timeout=30) as resp:
        html = resp.read().decode()

    m = re.search(r'RID = (\w+)', html)
    return m.group(1) if m else None


def _poll_blast(rid: str, max_wait: int = 180) -> str | None:
    params = urllib.parse.urlencode({
        "CMD":         "Get",
        "RID":         rid,
        "FORMAT_TYPE": "XML",
    })
    url = f"{BLAST_URL}?{params}"

    for _ in range(max_wait // 10):
        time.sleep(10)
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=30) as resp:
            content = resp.read().decode()

        if "Status=WAITING" in content:
            continue
        if "Status=FAILED" in content:
            return None
        if "Status=READY" in content or "<BlastOutput>" in content:
            return content

    return None


def _parse_blast_xml(xml_text: str) -> dict | None:
    from io import StringIO
    from Bio.Blast import NCBIXML

    try:
        blast_records = list(NCBIXML.parse(StringIO(xml_text)))
    except Exception:
        return None

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
        "chromosome":     chrom,
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
        "strand":         strand,
        "mismatches":     mismatch_count,
        "raw_mismatches": raw_mismatches,
        "query_length":   blast_record.query_length,
        "accession":      alignment.accession,
        "hit_title":      alignment.title,
    }


def _extract_chromosome(title: str) -> str:
    m = re.search(r'chromosome\s+(\d+|X|Y|MT|M)\b', title, re.IGNORECASE)
    if m:
        return m.group(1).upper()

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
