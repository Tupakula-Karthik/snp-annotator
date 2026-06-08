"""
dbsnp_lookup.py
───────────────
For each mismatch position from BLAST, queries NCBI Entrez/dbSNP to:
  1. Find rs IDs at that genomic coordinate
  2. Fetch variant details: gene, consequence, clinical significance

Uses Entrez eSearch + eSummary — no API key required (anonymous, 3 req/sec).
"""

import time
import re
from Bio import Entrez

# NCBI requires an email for anonymous Entrez access
Entrez.email = "snp-annotator-app@example.com"
RATE_LIMIT_DELAY = 0.4   # seconds between requests (≤3/sec anonymous)


def lookup_dbsnp(blast_result: dict) -> list[dict]:
    """
    Parameters
    ----------
    blast_result : dict from blast_aligner.run_blast()

    Returns
    -------
    list of variant dicts, one per mismatch:
    {
        position            : int
        ref_base            : str
        query_base          : str
        rs_id               : str | None      e.g. "rs1801133"
        gene                : str | None
        consequence         : str | None
        clinical_significance: str | None
    }
    """
    chrom      = blast_result.get("chromosome", "?")
    mismatches = blast_result.get("raw_mismatches", [])

    if not mismatches:
        return []

    variants = []
    for mm in mismatches:
        pos        = mm["position"]
        query_base = mm["query_base"]
        ref_base   = mm["ref_base"]

        ann = _query_single_position(chrom, pos)
        variants.append({
            "position":             pos,
            "ref_base":             ref_base,
            "query_base":           query_base,
            "rs_id":                ann.get("rs_id"),
            "gene":                 ann.get("gene"),
            "consequence":          ann.get("consequence"),
            "clinical_significance":ann.get("clinical_significance"),
        })
        time.sleep(RATE_LIMIT_DELAY)

    return variants


def _query_single_position(chrom: str, position: int) -> dict:
    """
    Search dbSNP for an rs ID at the given chromosome:position.
    Falls back gracefully at each step.
    """
    try:
        rs_id = _esearch_snp(chrom, position)
    except Exception:
        rs_id = None

    if not rs_id:
        return {}

    try:
        details = _esummary_snp(rs_id)
    except Exception:
        details = {}

    return {
        "rs_id": rs_id,
        **details,
    }


def _esearch_snp(chrom: str, position: int) -> str | None:
    """
    Use Entrez eSearch on the SNP database to find variants at a coordinate.
    Search term:  "chrN[CHR] AND pos[CHRPOS]"
    """
    # Normalise chromosome notation
    chrom_term = _normalise_chrom(chrom)
    term = f'{chrom_term}[CHR] AND {position}[CHRPOS]'

    handle  = Entrez.esearch(db="snp", term=term, retmax=1)
    record  = Entrez.read(handle)
    handle.close()
    time.sleep(RATE_LIMIT_DELAY)

    ids = record.get("IdList", [])
    if not ids:
        return None

    return f"rs{ids[0]}"


def _esummary_snp(rs_id: str) -> dict:
    """
    Fetch a DocSum for the given rs ID and extract key annotation fields.
    """
    numeric_id = rs_id.lstrip("rs")
    handle = Entrez.esummary(db="snp", id=numeric_id)
    records = Entrez.read(handle)
    handle.close()
    time.sleep(RATE_LIMIT_DELAY)

    # DocSum structure varies; try common paths
    doc = {}
    if isinstance(records, list) and records:
        doc = records[0]
    elif isinstance(records, dict):
        doc = records

    # DocumentSummarySet path
    if "DocumentSummarySet" in doc:
        summaries = doc["DocumentSummarySet"].get("DocumentSummary", [])
        if summaries:
            doc = summaries[0]

    return _parse_docsummary(doc)


def _parse_docsummary(doc: dict) -> dict:
    """Extract gene, consequence, and clinical significance from a DocSum."""
    result = {
        "gene":                 None,
        "consequence":          None,
        "clinical_significance":None,
    }

    # Gene name — various possible keys
    for key in ("GENES", "gene", "GENE_ID"):
        val = doc.get(key)
        if val:
            if isinstance(val, list) and val:
                gene_entry = val[0]
                result["gene"] = (
                    gene_entry.get("NAME") or
                    gene_entry.get("name") or
                    str(gene_entry)
                )[:30]
            elif isinstance(val, str) and val.strip():
                result["gene"] = val.strip()[:30]
            break

    # Functional consequence
    for key in ("FXN_CLASS", "fxn_class", "CONSEQUENCE"):
        val = doc.get(key)
        if val:
            result["consequence"] = _humanise_consequence(str(val))
            break

    # Clinical significance (ClinVar-linked)
    for key in ("CLINICAL_SIGNIFICANCE", "clinical_significance", "CLINSIG"):
        val = doc.get(key)
        if val and str(val).strip() not in ("", "0", "unknown"):
            result["clinical_significance"] = _humanise_clinsig(str(val))
            break

    return result


def _normalise_chrom(chrom: str) -> str:
    """Convert various chromosome notations to the integer/letter form dbSNP expects."""
    chrom = str(chrom).upper().replace("CHR", "").strip()
    mapping = {"MT": "26", "M": "26", "X": "23", "Y": "24"}
    return mapping.get(chrom, chrom)


def _humanise_consequence(raw: str) -> str:
    """Map dbSNP FXN_CLASS codes to human-readable strings."""
    mapping = {
        "missense":         "Missense variant",
        "nonsense":         "Nonsense (stop-gain)",
        "synonymous":       "Synonymous (silent)",
        "frameshift":       "Frameshift",
        "intron":           "Intronic",
        "splice":           "Splice site",
        "utr":              "UTR variant",
        "downstream":       "Downstream gene variant",
        "upstream":         "Upstream gene variant",
        "intergenic":       "Intergenic",
        "cds":              "Coding sequence variant",
        "3prime_utr":       "3′ UTR variant",
        "5prime_utr":       "5′ UTR variant",
    }
    raw_lower = raw.lower()
    for key, label in mapping.items():
        if key in raw_lower:
            return label
    return raw.replace("_", " ").title()


def _humanise_clinsig(raw: str) -> str:
    """Normalise ClinVar clinical significance strings."""
    raw_lower = raw.lower()
    if "pathogenic" in raw_lower and "likely" not in raw_lower:
        return "Pathogenic"
    if "likely_pathogenic" in raw_lower or "likely pathogenic" in raw_lower:
        return "Likely pathogenic"
    if "likely_benign" in raw_lower or "likely benign" in raw_lower:
        return "Likely benign"
    if "benign" in raw_lower:
        return "Benign"
    if "uncertain" in raw_lower or "vus" in raw_lower:
        return "Uncertain significance"
    if "risk" in raw_lower:
        return "Risk factor"
    return raw.replace("_", " ").title()
