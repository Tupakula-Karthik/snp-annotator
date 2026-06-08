"""
ab1_parser.py
─────────────
Reads an ABI .ab1 Sanger sequencing file using BioPython and returns:
  - seq_record   : SeqRecord with the called base sequence + Phred quality scores
  - quality      : list[int] of Phred scores (trimmed to good region)
  - raw_channels : dict of raw trace channel arrays (A/C/G/T) for chromatogram
"""

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np


def parse_ab1(file_path: Path) -> tuple[SeqRecord, list[int], dict]:
    """
    Parameters
    ----------
    file_path : Path  –  path to .ab1 file

    Returns
    -------
    seq_record     : BioPython SeqRecord (sequence + quality in .letter_annotations)
    quality_scores : list of Phred scores (same length as sequence)
    raw_channels   : dict {"A":[], "C":[], "G":[], "T":[]} of raw trace intensities
    """
    record = SeqIO.read(str(file_path), "abi")

    # Quality scores
    quality = record.letter_annotations.get("phred_quality", [])
    if not quality:
        # Fall back: try PCON1 tag
        quality = list(record.annotations.get("abif_raw", {}).get("PCON1", []))
    if not quality:
        quality = [30] * len(record.seq)   # placeholder if unavailable

    # Raw trace channels (used for chromatogram rendering)
    raw = record.annotations.get("abif_raw", {})
    channels = {}
    channel_keys = {
        "A": ["DATA9",  "DATA1"],
        "C": ["DATA10", "DATA2"],
        "G": ["DATA11", "DATA3"],
        "T": ["DATA12", "DATA4"],
    }
    for base, keys in channel_keys.items():
        for k in keys:
            if k in raw:
                channels[base] = list(raw[k])
                break

    # Trim low-quality ends (Phred < 20)
    trimmed_record, trimmed_quality = _trim_quality(record, quality, cutoff=20)

    return trimmed_record, trimmed_quality, channels


def _trim_quality(record: SeqRecord, quality: list, cutoff: int = 20):
    """
    Remove low-quality bases from the 5' and 3' ends using a sliding window.
    Returns trimmed SeqRecord and trimmed quality list.
    """
    if len(quality) < 20:
        return record, quality

    window = 10
    start = 0
    for i in range(len(quality) - window):
        if np.mean(quality[i:i+window]) >= cutoff:
            start = i
            break

    end = len(quality)
    for i in range(len(quality) - window, 0, -1):
        if np.mean(quality[i:i+window]) >= cutoff:
            end = i + window
            break

    trimmed_seq    = record.seq[start:end]
    trimmed_qual   = quality[start:end]

    new_record = record[start:end]
    new_record.letter_annotations["phred_quality"] = trimmed_qual

    return new_record, trimmed_qual
