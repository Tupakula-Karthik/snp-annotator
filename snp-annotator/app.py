import streamlit as st
import pandas as pd
import time
from pathlib import Path

from pipeline.ab1_parser import parse_ab1
from pipeline.blast_aligner import run_blast
from pipeline.dbsnp_lookup import lookup_dbsnp
from pipeline.utils import render_chromatogram

st.set_page_config(
    page_title="SNP Annotator",
    page_icon="🧬",
    layout="wide",
)

# ── Styles ──────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    .block-container { padding-top: 2rem; }
    .stProgress > div > div { background-color: #4CAF50; }
    .step-card {
        border-left: 4px solid #4CAF50;
        padding: 0.6rem 1rem;
        margin: 0.4rem 0;
        background: #f8f9fa;
        border-radius: 0 6px 6px 0;
    }
    .step-card.pending { border-color: #ccc; }
    .step-card.running { border-color: #2196F3; background: #e3f2fd; }
    .step-card.done    { border-color: #4CAF50; background: #e8f5e9; }
    .step-card.error   { border-color: #f44336; background: #ffebee; }
    .snp-badge {
        display: inline-block;
        padding: 2px 8px;
        border-radius: 12px;
        font-size: 0.8rem;
        font-weight: 600;
    }
    .badge-path { background: #ffcdd2; color: #c62828; }
    .badge-likely { background: #ffe0b2; color: #e65100; }
    .badge-benign { background: #c8e6c9; color: #1b5e20; }
    .badge-unknown { background: #e0e0e0; color: #424242; }
</style>
""", unsafe_allow_html=True)

# ── Header ───────────────────────────────────────────────────────────────────
col1, col2 = st.columns([3, 1])
with col1:
    st.title("🧬 Sanger SNP Annotator")
    st.caption("Upload an .ab1 file → align to genome → find & annotate SNP variants")
with col2:
    st.markdown("<br>", unsafe_allow_html=True)
    st.info("Anonymous NCBI access\n\n~3 req/sec limit applies", icon="ℹ️")

st.divider()

# ── File Upload ───────────────────────────────────────────────────────────────
uploaded = st.file_uploader(
    "Upload your .ab1 Sanger sequencing file",
    type=["ab1"],
    help="ABI sequencing files typically from Applied Biosystems instruments",
)

if not uploaded:
    st.markdown("""
    ### How it works
    1. **Parse** — BioPython extracts the base sequence and Phred quality scores from your .ab1 file
    2. **Align** — NCBI BLAST maps your sequence to the human reference genome (GRCh38) and returns chromosomal coordinates
    3. **Annotate** — NCBI Entrez queries dbSNP at every position that differs from the reference, returning rs IDs and clinical significance
    """)
    st.stop()

# ── Pipeline ──────────────────────────────────────────────────────────────────
run_btn = st.button("▶  Run Analysis", type="primary", use_container_width=True)

if not run_btn:
    st.info("Click **Run Analysis** to start the pipeline.", icon="👆")
    st.stop()

progress_bar = st.progress(0)
step_slots = {}
for key, label in [("parse","Step 1 — Parse .ab1"), ("blast","Step 2 — BLAST alignment"), ("dbsnp","Step 3 — dbSNP lookup"), ("results","Step 4 — Render results")]:
    step_slots[key] = st.empty()
    step_slots[key].markdown(f'<div class="step-card pending">⏳ {label}</div>', unsafe_allow_html=True)

log_slot = st.empty()

def update(key, label, status, msg=""):
    icons = {"running":"🔄","done":"✅","error":"❌"}
    step_slots[key].markdown(
        f'<div class="step-card {status}">{icons.get(status,"⏳")} {label}{"  —  "+msg if msg else ""}</div>',
        unsafe_allow_html=True
    )

# Step 1 — Parse
update("parse","Step 1 — Parse .ab1","running")
progress_bar.progress(5)
try:
    tmp_path = Path("/tmp") / uploaded.name
    tmp_path.write_bytes(uploaded.getvalue())
    seq_record, quality_scores, raw_bases = parse_ab1(tmp_path)
    update("parse","Step 1 — Parse .ab1","done", f"{len(seq_record.seq)} bp | avg quality {sum(quality_scores)//len(quality_scores)}")
    progress_bar.progress(25)
except Exception as e:
    update("parse","Step 1 — Parse .ab1","error", str(e))
    st.error(f"Failed to parse .ab1 file: {e}")
    st.stop()

# Show chromatogram preview
with st.expander("📊 Chromatogram preview (first 200 bp)", expanded=False):
    fig = render_chromatogram(seq_record, raw_bases, limit=200)
    if fig:
        st.pyplot(fig)
    else:
        st.caption("Raw trace data not available in this file.")

# Step 2 — BLAST
update("blast","Step 2 — BLAST alignment","running")
progress_bar.progress(30)
log_slot.caption("⏳ Submitting to NCBI BLAST… this typically takes 30–90 seconds on anonymous access.")
try:
    blast_result = run_blast(str(seq_record.seq))
    if not blast_result:
        raise ValueError("No BLAST hits returned. Sequence may be too short or non-human.")
    chrom  = blast_result["chromosome"]
    hsp    = blast_result["hsp"]
    strand = blast_result["strand"]
    mismatches = blast_result["mismatches"]
    update("blast","Step 2 — BLAST alignment","done",
           f"chr{chrom}:{hsp['sbjct_start']:,}–{hsp['sbjct_end']:,}  |  {mismatches} mismatch(es)  |  strand {'+'if strand==1 else '-'}")
    progress_bar.progress(60)
    log_slot.empty()
except Exception as e:
    update("blast","Step 2 — BLAST alignment","error", str(e))
    st.error(f"BLAST alignment failed: {e}")
    st.stop()

# Step 3 — dbSNP
update("dbsnp","Step 3 — dbSNP lookup","running")
progress_bar.progress(65)
log_slot.caption(f"⏳ Querying dbSNP for {len(mismatches) if isinstance(mismatches, list) else mismatches} candidate variant position(s)…")
try:
    variants = lookup_dbsnp(blast_result)
    update("dbsnp","Step 3 — dbSNP lookup","done", f"{len(variants)} variant(s) annotated")
    progress_bar.progress(90)
    log_slot.empty()
except Exception as e:
    update("dbsnp","Step 3 — dbSNP lookup","error", str(e))
    st.warning(f"dbSNP lookup encountered an issue: {e}. Showing unannotated mismatches.")
    variants = blast_result.get("raw_mismatches", [])

# Step 4 — Results
update("results","Step 4 — Render results","running")
progress_bar.progress(95)

def sig_badge(sig):
    if not sig or sig in ("unknown", "not_provided", ""):
        return '<span class="snp-badge badge-unknown">Unknown</span>'
    s = sig.lower()
    if "pathogenic" in s and "likely" not in s:
        return f'<span class="snp-badge badge-path">{sig}</span>'
    if "likely" in s:
        return f'<span class="snp-badge badge-likely">{sig}</span>'
    if "benign" in s:
        return f'<span class="snp-badge badge-benign">{sig}</span>'
    return f'<span class="snp-badge badge-unknown">{sig}</span>'

progress_bar.progress(100)
update("results","Step 4 — Render results","done")

st.divider()
st.subheader(f"🔬 Results — chr{chrom}:{hsp['sbjct_start']:,}–{hsp['sbjct_end']:,}")

if not variants:
    st.success("✅ No mismatches found relative to the reference genome. Your sequence perfectly matches GRCh38.")
else:
    # Summary metrics
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Total variants", len(variants))
    known = sum(1 for v in variants if v.get("rs_id") and v["rs_id"] != "—")
    m2.metric("Known (rs ID)", known)
    path  = sum(1 for v in variants if "pathogenic" in (v.get("clinical_significance") or "").lower())
    m3.metric("Pathogenic flags", path)
    m4.metric("Alignment region", f"chr{chrom}")

    st.markdown("<br>", unsafe_allow_html=True)

    # Build display table
    rows = []
    for v in variants:
        sig   = v.get("clinical_significance") or "unknown"
        rs    = v.get("rs_id") or "—"
        rslink = f'<a href="https://www.ncbi.nlm.nih.gov/snp/{rs}" target="_blank">{rs}</a>' if rs != "—" else "—"
        rows.append({
            "Position (GRCh38)": f"chr{chrom}:{v.get('position', '?'):,}" if isinstance(v.get('position'), int) else f"chr{chrom}:{v.get('position','?')}",
            "Ref": v.get("ref_base", "?"),
            "Your base": v.get("query_base", "?"),
            "rs ID": rslink,
            "Gene": v.get("gene", "—"),
            "Consequence": v.get("consequence", "—"),
            "Clinical significance": sig_badge(sig),
        })

    df_html = pd.DataFrame(rows).to_html(escape=False, index=False, classes="dataframe")
    st.markdown(df_html, unsafe_allow_html=True)

    # Download CSV
    csv_rows = []
    for v in variants:
        csv_rows.append({
            "chromosome": f"chr{chrom}",
            "position": v.get("position",""),
            "ref_base": v.get("ref_base",""),
            "query_base": v.get("query_base",""),
            "rs_id": v.get("rs_id",""),
            "gene": v.get("gene",""),
            "consequence": v.get("consequence",""),
            "clinical_significance": v.get("clinical_significance",""),
        })
    csv_data = pd.DataFrame(csv_rows).to_csv(index=False)
    st.download_button(
        "⬇️  Download results as CSV",
        data=csv_data,
        file_name=f"snp_results_{uploaded.name.replace('.ab1','')}.csv",
        mime="text/csv",
    )

st.divider()
st.caption("Data sources: NCBI BLAST (GRCh38 alignment) · NCBI dbSNP (variant annotation) · ClinVar (clinical significance) — Anonymous access, rate limits apply.")
