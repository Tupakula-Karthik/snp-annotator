# 🧬 Sanger SNP Annotator

Upload an `.ab1` Sanger sequencing file and automatically find SNP variations with genomic annotations.

## What it does

| Step | Tool | Output |
|------|------|--------|
| Parse `.ab1` | BioPython | Base sequence + Phred quality scores |
| Align to genome | NCBI BLAST (GRCh38) | Chromosome + genomic coordinates |
| Find SNPs | Mismatch detection | Position, ref base, query base |
| Annotate | NCBI Entrez / dbSNP | rs IDs, gene, consequence, clinical significance |

---

## Deploy to Streamlit Cloud (free, shareable URL)

### Prerequisites
- GitHub account (free)
- Streamlit Community Cloud account (free) → [share.streamlit.io](https://share.streamlit.io)

### Steps

**1. Push to GitHub**
```bash
# Create a new repo on github.com, then:
git init
git add .
git commit -m "Initial SNP annotator"
git remote add origin https://github.com/YOUR_USERNAME/snp-annotator.git
git push -u origin main
```

**2. Deploy on Streamlit Cloud**
1. Go to [share.streamlit.io](https://share.streamlit.io) → **New app**
2. Connect your GitHub account
3. Select your repo → branch `main` → main file `app.py`
4. Click **Deploy**

That's it. In ~2 minutes you get a URL like:
```
https://your-username-snp-annotator-app-xxxxxx.streamlit.app
```
Share this link with anyone — no installation needed on their end.

---

## Optional: Add your NCBI API key (recommended)

Anonymous NCBI access is limited to ~3 requests/second and can be slow.
A free API key raises the limit to 10/second.

1. Get a key at [ncbi.nlm.nih.gov/account](https://www.ncbi.nlm.nih.gov/account/)
2. In Streamlit Cloud: **App settings → Secrets**
3. Add:
```toml
NCBI_API_KEY = "your_key_here"
```
4. The app will automatically detect and use it.

---

## Run locally (optional)

```bash
pip install -r requirements.txt
streamlit run app.py
```

---

## File structure

```
snp-annotator/
├── app.py                  # Main Streamlit interface
├── requirements.txt        # Python dependencies
├── .streamlit/
│   └── config.toml         # Theme + server settings
└── pipeline/
    ├── __init__.py
    ├── ab1_parser.py        # BioPython .ab1 parsing
    ├── blast_aligner.py     # NCBI BLAST alignment
    ├── dbsnp_lookup.py      # Entrez / dbSNP annotation
    └── utils.py             # Chromatogram renderer
```

---

## Notes on BLAST timing

NCBI BLAST on anonymous access typically takes **30–90 seconds** per sequence.
The app shows a live progress indicator while waiting.

## Data sources

- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) — alignment against GRCh38
- [NCBI dbSNP](https://www.ncbi.nlm.nih.gov/snp/) — variant rs IDs
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) — clinical significance (via dbSNP links)
