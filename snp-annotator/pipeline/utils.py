"""
utils.py
────────
Helper utilities:
  - render_chromatogram : draws the raw trace channels as a matplotlib figure
"""

import matplotlib
matplotlib.use("Agg")   # non-interactive backend for Streamlit
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from Bio.SeqRecord import SeqRecord


BASE_COLORS = {
    "A": "#2ecc71",   # green
    "C": "#3498db",   # blue
    "G": "#e67e22",   # orange
    "T": "#e74c3c",   # red
}


def render_chromatogram(seq_record: SeqRecord, raw_channels: dict, limit: int = 200):
    """
    Draw the first `limit` called bases with their underlying trace signals.

    Parameters
    ----------
    seq_record   : parsed BioPython SeqRecord
    raw_channels : dict {"A":list, "C":list, "G":list, "T":list}
    limit        : number of bases / trace points to show

    Returns
    -------
    matplotlib Figure or None if trace data unavailable
    """
    if not raw_channels or len(raw_channels) < 4:
        return None

    # Normalise channels
    channels = {}
    for base, data in raw_channels.items():
        arr = np.array(data[:limit * 10], dtype=float)  # rough estimate
        if arr.max() > 0:
            arr = arr / arr.max()
        channels[base] = arr

    seq_str = str(seq_record.seq)[:limit]
    quality = seq_record.letter_annotations.get("phred_quality", [])[:limit]

    fig, (ax_trace, ax_qual) = plt.subplots(
        2, 1,
        figsize=(14, 4),
        gridspec_kw={"height_ratios": [3, 1]},
    )
    fig.patch.set_facecolor("#0e1117")
    for ax in (ax_trace, ax_qual):
        ax.set_facecolor("#0e1117")
        for spine in ax.spines.values():
            spine.set_color("#444")
        ax.tick_params(colors="#aaa")

    # ── Trace plot ──────────────────────────────────────────
    # Determine trace length to plot
    trace_len = min(len(v) for v in channels.values())
    xs = np.linspace(0, len(seq_str), trace_len)

    for base in ("A", "C", "G", "T"):
        if base in channels:
            ys = channels[base][:trace_len]
            ax_trace.plot(xs, ys, color=BASE_COLORS[base], linewidth=0.8, alpha=0.85, label=base)

    # Overlay base calls
    for i, base in enumerate(seq_str):
        x_pos = i + 0.5
        ax_trace.text(
            x_pos, 1.02, base,
            ha="center", va="bottom",
            fontsize=6.5,
            color=BASE_COLORS.get(base.upper(), "#fff"),
            fontfamily="monospace",
            transform=ax_trace.get_xaxis_transform(),
        )

    ax_trace.set_xlim(0, len(seq_str))
    ax_trace.set_ylim(-0.05, 1.15)
    ax_trace.set_ylabel("Intensity", color="#aaa", fontsize=9)
    ax_trace.set_title("Chromatogram (first 200 bp)", color="#ddd", fontsize=10, pad=12)

    legend_patches = [mpatches.Patch(color=BASE_COLORS[b], label=b) for b in "ACGT"]
    ax_trace.legend(handles=legend_patches, loc="upper right", fontsize=8,
                    facecolor="#222", labelcolor="#ddd", edgecolor="#444")

    # ── Quality plot ─────────────────────────────────────────
    if quality:
        q_arr = np.array(quality)
        colors = np.where(q_arr >= 30, "#2ecc71", np.where(q_arr >= 20, "#f39c12", "#e74c3c"))
        ax_qual.bar(np.arange(len(q_arr)) + 0.5, q_arr, color=colors, width=0.9)
        ax_qual.axhline(20, color="#f39c12", linewidth=0.7, linestyle="--", alpha=0.6)
        ax_qual.axhline(30, color="#2ecc71", linewidth=0.7, linestyle="--", alpha=0.6)
        ax_qual.set_xlim(0, len(seq_str))
        ax_qual.set_ylabel("Phred Q", color="#aaa", fontsize=9)
        ax_qual.set_xlabel("Position (bp)", color="#aaa", fontsize=9)

    plt.tight_layout(pad=0.5)
    return fig
