"""
dbsnp_lookup.py
"""

import time
import re
import json
import urllib.request
from Bio import Entrez

Entrez.email = "snp-annotator-app@example.com"
RATE_LIMIT_DELAY = 0.4


def lookup_dbsnp(blast_result: dict) -> list[dict]:
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
            "position":              pos,
            "ref_base":              ref_base,
            "query_base":            query_base,
            "rs_id":                 ann.get("rs_id"),
            "gene":                  ann.get("gene"),
            "consequence":           ann.get("consequence"),
            "clinical_significance": ann.get("clinical_significance"),
        })
        time.sleep(RATE_LIMIT_DELAY)

    return variants


def _query_single_position(chrom: str, position: int) -> dict:
    try:
        rs_id = _esearch_snp(chrom, position)
    except Exception:
        rs_id = None

    if not rs_id:
        return {}

    try:
        details = _fetch_snp_via_api(rs_id)
    except Exception:
        details = {}

    return {"rs_id": rs_id, **details}


def _esearch_snp(chrom: str, position: int) -> str | None:
    chrom_term = _normalise_chrom(chrom)
    term = f'{chrom_term}[CHR] AND {position}[CHRPOS]'

    handle = Entrez.esearch(db="snp", term=term, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    time.sleep(RATE_LIMIT_DELAY)

    ids = record.get("IdList", [])
    if not ids:
        return None

    return f"rs{ids[0]}"


def _fetch_snp_via_api(rs_id: str) -> dict:
    """Use NCBI Variation Services REST API — more reliable than eSummary."""
    numeric_id = rs_id.lstrip("rs")
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{numeric_id}"

    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
        return _parse_variation_api(data)
    except Exception:
        return _esummary_snp(rs_id)


def _parse_variation_api(data: dict) -> dict:
    result = {
        "gene":                  None,
        "consequence":           None,
        "clinical_significance": None,
    }

    try:
        allele_annotations = (
            data.get("primary_snapshot_data", {})
                .get("allele_annotations", [])
        )
        genes = set()
        consequences = set()

        for allele in allele_annotations:
            for ann in allele.get("assembly_annotation", []):
                for gene in ann.get("genes", []):
                    name = gene.get("name") or gene.get("locus")
                    if name:
                        genes.add(name)
                    for rna in gene.get("rnas", []):
                        for c in rna.get("sequence_ontology", []):
                            label = c.get("name", "")
                            if label:
                                consequences.add(_humanise_consequence(label))

        if genes:
            result["gene"] = ", ".join(sorted(genes)[:2])
        if consequences:
            result["consequence"] = ", ".join(sorted(consequences)[:2])

    except Exception:
        pass

    try:
        clin_sigs = set()
        support = (
            data.get("primary_snapshot_data", {})
                .get("support", [])
        )
        for s in support:
            for interp in s.get("clinical_significances", []):
                sig = interp.get("clinical_significance_type", "")
                if sig:
                    clin_sigs.add(_humanise_clinsig(sig))

        if clin_sigs:
            result["clinical_significance"] = ", ".join(sorted(clin_sigs))
    except Exception:
        pass

    return result


def _esummary_snp(rs_id: str) -> dict:
    numeric_id = rs_id.lstrip("rs")
    handle = Entrez.esummary(db="snp", id=numeric_id)
    records = Entrez.read(handle)
    handle.close()
    time.sleep(RATE_LIMIT_DELAY)

    doc = {}
    if isinstance(records, list) and records:
        doc = records[0]
    elif isinstance(records, dict):
        doc = records

    if "DocumentSummarySet" in doc:
        summaries = doc["DocumentSummarySet"].get("DocumentSummary", [])
        if summaries:
            doc = summaries[0]

    return _parse_docsummary(doc)


def _parse_docsummary(doc: dict) -> dict:
    result = {
        "gene":                  None,
        "consequence":           None,
        "clinical_significance": None,
    }

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

    for key in ("FXN_CLASS", "fxn_class", "CONSEQUENCE"):
        val = doc.get(key)
        if val:
            result["consequence"] = _humanise_consequence(str(val))
            break

    for key in ("CLINICAL_SIGNIFICANCE", "clinical_significance", "CLINSIG"):
        val = doc.get(key)
        if val and str(val).strip() not in ("", "0", "unknown"):
            result["clinical_significance"] = _humanise_clinsig(str(val))
            break

    return result


def _normalise_chrom(chrom: str) -> str:
    chrom = str(chrom).upper().replace("CHR", "").strip()
    mapping = {"MT": "26", "M": "26", "X": "23", "Y": "24"}
    return mapping.get(chrom, chrom)


def _humanise_consequence(raw: str) -> str:
    mapping = {
        "missense":          "Missense variant",
        "nonsense":          "Nonsense (stop-gain)",
        "synonymous":        "Synonymous (silent)",
        "frameshift":        "Frameshift",
        "intron":            "Intronic",
        "splice":            "Splice site",
        "utr":               "UTR variant",
        "downstream":        "Downstream gene variant",
        "upstream":          "Upstream gene variant",
        "intergenic":        "Intergenic",
        "cds":               "Coding sequence variant",
        "3prime_utr":        "3′ UTR variant",
        "5prime_utr":        "5′ UTR variant",
        "stop_gained":       "Nonsense (stop-gain)",
        "stop_lost":         "Stop lost",
        "start_lost":        "Start lost",
        "nc_transcript":     "Non-coding transcript variant",
    }
    raw_lower = raw.lower()
    for key, label in mapping.items():
        if key in raw_lower:
            return label
    return raw.replace("_", " ").title()


def _humanise_clinsig(raw: str) -> str:
    raw_lower = raw.lower()
    if "likely_pathogenic" in raw_lower or "likely pathogenic" in raw_lower:
        return "Likely pathogenic"
    if "pathogenic" in raw_lower:
        return "Pathogenic"
    if "likely_benign" in raw_lower or "likely benign" in raw_lower:
        return "Likely benign"
    if "benign" in raw_lower:
        return "Benign"
    if "uncertain" in raw_lower or "vus" in raw_lower:
        return "Uncertain significance"
    if "risk" in raw_lower:
        return "Risk factor"
    return raw.replace("_", " ").title()
