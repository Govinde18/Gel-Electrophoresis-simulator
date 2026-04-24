"""
Gel Electrophoresis Simulator v8.0
"""

# ── Standard library ──────────────────────────────────────────────────────────
import os
import tempfile
import threading
import numpy as np
from datetime import datetime

import tkinter as tk
from tkinter import filedialog as _filedialog

ROOT = os.path.dirname(os.path.abspath(__file__))

# ── Matplotlib (headless backend for Kivy embedding) ──────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ── Biopython ─────────────────────────────────────────────────────────────────
from Bio.Restriction import AllEnzymes, Analysis
from Bio.Seq import Seq
try:
    from Bio import SeqIO as _SeqIO
    _HAS_SEQIO = True
except ImportError:
    _HAS_SEQIO = False

# ── Kivy ──────────────────────────────────────────────────────────────────────
from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.gridlayout import GridLayout
from kivy.uix.scrollview import ScrollView
from kivy.uix.button import Button
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.spinner import Spinner
from kivy.uix.checkbox import CheckBox
from kivy.uix.slider import Slider
from kivy.uix.image import Image
from kivy.uix.popup import Popup
from kivy.core.window import Window
from kivy.metrics import dp
from kivy.clock import Clock
from kivy.graphics import Color, Rectangle, RoundedRectangle

# ══════════════════════════════════════════════════════════════════════════════
#  THEME
# ══════════════════════════════════════════════════════════════════════════════
_DARK     = (0.11, 0.11, 0.14, 1)
_PANEL    = (0.17, 0.17, 0.22, 1)
_ACCENT   = (0.30, 0.68, 0.98, 1)   # blue
_GREEN    = (0.30, 0.82, 0.52, 1)   # run / success
_DANGER   = (0.93, 0.35, 0.35, 1)   # error / remove
_BTN      = (0.24, 0.24, 0.32, 1)
_TEXT     = (0.94, 0.94, 0.96, 1)
_HINT     = (0.56, 0.57, 0.66, 1)

Window.clearcolor = _DARK

# ══════════════════════════════════════════════════════════════════════════════
#  CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════
PRESET_LADDERS = {
    "1 kb Plus – NEB":         [10000,8000,6000,5000,4000,3000,2500,2000,1500,1000,750,500,250,100],
    "1 kb – NEB":              [10000,8000,6000,5000,4000,3000,2500,2000,1500,1000,750,500],
    "100 bp – NEB":            [1000,900,800,700,600,500,400,300,200,100],
    "GeneRuler 1 kb Plus":     [20000,10000,7000,5000,4000,3000,2000,1500,1000,700,500,400,300,200,75],
    "λ-EcoRI/BglII (Takara)":  [22010,19329,13286,9688,7743,6223,4254,3472,2392,1882,1489,925,651,421,415,74],
    "Custom (enter below)":    [],
}

GEL_OPTS    = ["0.5%","0.8%","1.0%","1.2%","1.5%","2.0%","3.0%"]
BUF_OPTS    = ["TAE","TBE"]

# Buffer biophysical factors (speed, resolution)
BUF_FACTORS = {
    "TAE": {"speed": 1.00, "resolution": 1.00},
    "TBE": {"speed": 0.82, "resolution": 1.15},   # TBE: slower but sharper bands
}


# Enzyme names are taken verbatim from Biopython's REBASE-derived AllEnzymes
# so Biopython's Analysis engine applies the correct recognition+cut site for each.
# Recognition sequences and canonical cut sites shown for reference:
#   EcoRI   G^AATTC       BamHI   G^GATCC     HindIII  A^AGCTT
#   NcoI    C^CATGG       NdeI    CA^TATG      XhoI     C^TCGAG
#   SalI    G^TCGAC       KpnI    GGTAC^C      SacI     GAGCT^C
#   XbaI    T^CTAGA       SphI    GCATG^C      PstI     CTGCA^G
#   SmaI    CCC^GGG       MluI    A^CGCGT      NheI     G^CTAGC
#   ClaI    AT^CGAT       EcoRV   GAT^ATC      NotI     GC^GGCCGC
#   BglII   A^GATCT       AgeI    A^CCGGT      NsiI     ATGCA^T
#   AvrII   C^CTAGG       SpeI    A^CTAGT      NarI     GG^CGCC
#   StuI    AGG^CCT       SacII   CCGC^GG      MfeI     C^AATTG
#   AflII   C^TTAAG       ApoI    R^AATTY      AscI     GG^CGCGCC
#   PacI    TTAAT^TAA     SwaI    ATTT^AAAT    FseI     GGCCGG^CC
#   AsiSI   GCGAT^CGC     SgrAI   CR^CCGGYG    BstBI    TT^CGAA
#   BspHI   T^CATGA       BspEI   T^CCGGA      BsrGI    T^GTACA
#   BstEII  G^GTNACC      BstXI   CCANN^NNNTGG MscI     TGG^CCA
#   NgoMIV  G^CCGGC       PmeI    GTTT^AAAC    SfiI     GGCCN^NNNNGGCC
#   SrfI    GCCC^GGGC     StuI    AGG^CCT      Tth111I  GACN^NNGTC
#   XmaI    C^CCGGG       BsiWI   C^GTACG      BlpI     GC^TNAGC
#   BsaBI   GATN^NNNATC   BsaBI   (above)      BclI     T^GATCA
#   BseYI   CCCAGC(-5/-1) BsiHKAI GWGCW^C      BsmBI    CGTCTC(1/5)
#   BsrFI   R^CCGGY       BssHII  G^CGCGC      BstAPI   GCANNNN^NTGC
#   DraI    TTT^AAA       EagI    C^GGCCG       EcoNI    CCTNN^NNNAGG
#   HpaI    GTT^AAC       MluCI   ^AATT         PciI     A^CATGT
#   RsrII   CG^GWCCG      SbfI    CCTGCA^GG    SexAI    A^CCWGGTY
#   SmlI    C^TYRAG       SnaBI   TAC^GTA       SrfI     GCCC^GGGC
#   TaiI    ACGT           XcmI   CCANN^NNNNNNNNNTGG
COMMON_ENZYMES = sorted([
    # ── Workhorse 6-cutters ────────────────────────────────────────────────
    "EcoRI","BamHI","HindIII","NcoI","NdeI","XhoI","SalI",
    "KpnI","SacI","XbaI","SphI","PstI","SmaI","MluI","NheI",
    "ClaI","EcoRV","NotI","BglII","AgeI","NsiI","AvrII",
    # ── Common cloning / expression vector enzymes ─────────────────────────
    "SpeI","NarI","StuI","SacII","MfeI","AflII","AscI",
    "PacI","SwaI","FseI","AsiSI","SgrAI","BstBI","BspHI",
    "BspEI","BsrGI","BstEII","MscI","NgoMIV","PmeI","SfiI",
    "SrfI","XmaI","BsiWI","BlpI","BclI","BssHII","EagI",
    "HpaI","PciI","SbfI","SnaBI","DraI","EcoNI","RsrII",
    # ── Commonly used 4-cutters (higher frequency, useful for mapping) ─────
    "TaiI","MboI","Sau3AI","MseI","NlaIII","AluI","HaeIII",
    "RsaI","TaqI","MspI","HpaII","CfoI","HhaI","HinfI","BfaI",
    "MaeII","BstNI","DdeI","MboII","NlaIV","Tsp45I",
    # ── Rare-cutters (8-bp recognition) ────────────────────────────────────
    "SgfI","PacI","SwaI","AscI","FseI","PmeI","SrfI",
    # ── Additional common lab enzymes ──────────────────────────────────────
    "ApaI","ApaLI","BclI","BglI","BsaBI","BseYI","BsiHKAI",
    "BsmBI","BsrFI","BstAPI","BstXI","BstZ17I","Cfr10I","CpoI",
    "DrdI","DraIII","EarI","EcoT22I","Eco47III","EcoO109I",
    "Esp3I","HincII","HindII","KasI","MreI","MroI","MunI",
    "NaeI","NspI","NspBII","PflMI","PflFI","PmlI","PscI",
    "PspOMI","PsrI","PvuI","PvuII","SalI","ScaI","SciI",
    "SexAI","SgfI","SmlI","StuI","SwaI",
    "TfiI","Tth111I","Van91I","XcmI","XmnI","ZraI",
])

# ══════════════════════════════════════════════════════════════════════════════
#  FILE PARSING
# ══════════════════════════════════════════════════════════════════════════════
_DNA_CHARS = set("ACGTN")

def _clean_seq(raw: str) -> str:
    return "".join(c for c in raw.upper() if c in _DNA_CHARS)


def _parse_fasta(filepath: str) -> dict:
    """Return {id: sequence} for every record in a FASTA file."""
    if _HAS_SEQIO:
        try:
            return {r.id: str(r.seq) for r in _SeqIO.parse(filepath, "fasta")}
        except Exception:
            pass
    # Manual fallback
    records, name, buf = {}, None, []
    try:
        with open(filepath, "r", errors="ignore") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    if name:
                        records[name] = "".join(buf)
                    name = line[1:].split()[0]
                    buf = []
                elif name:
                    buf.append(line)
            if name:
                records[name] = "".join(buf)
    except Exception:
        pass
    return records


def _parse_genbank(filepath: str) -> dict:
    if _HAS_SEQIO:
        try:
            return {r.id: str(r.seq) for r in _SeqIO.parse(filepath, "genbank")}
        except Exception:
            pass
    return {}


def _parse_snapgene_dna(filepath: str) -> dict:
    """
    Parse a SnapGene .dna binary file.
    Format: repeated segments [1-byte type][4-byte big-endian length][data].
    Segment type 0x00 = DNA; first data byte = topology flags; rest = sequence.
    """
    try:
        with open(filepath, "rb") as fh:
            data = fh.read()
        idx = 0
        while idx + 5 <= len(data):
            seg_type = data[idx]
            seg_len  = int.from_bytes(data[idx+1:idx+5], "big")
            end      = idx + 5 + seg_len
            if end > len(data):
                break
            if seg_type == 0x00 and seg_len > 1:
                seq = _clean_seq(data[idx+6:end].decode("ascii", errors="ignore"))
                if len(seq) > 10:
                    name = os.path.splitext(os.path.basename(filepath))[0]
                    return {name: seq}
            idx = end
    except Exception:
        pass
    # Fallback: treat file as plain text
    try:
        with open(filepath, "r", errors="ignore") as fh:
            seq = _clean_seq(fh.read())
        if len(seq) > 10:
            name = os.path.splitext(os.path.basename(filepath))[0]
            return {name: seq}
    except Exception:
        pass
    return {}


def parse_sequence_file(filepath: str) -> dict:
    """Dispatch to the right parser; return {name: raw_sequence}."""
    ext = os.path.splitext(filepath)[1].lower()
    if ext in (".fa", ".fasta", ".fas"):
        result = _parse_fasta(filepath)
        if result:
            return result
    elif ext in (".gb", ".gbk", ".genbank"):
        result = _parse_genbank(filepath)
        if result:
            return result
    elif ext == ".dna":
        return _parse_snapgene_dna(filepath)
    # Generic: try FASTA, then plain sequence
    result = _parse_fasta(filepath)
    if result:
        return result
    try:
        with open(filepath, "r", errors="ignore") as fh:
            seq = _clean_seq(fh.read())
        if len(seq) > 10:
            name = os.path.splitext(os.path.basename(filepath))[0]
            return {name: seq}
    except Exception:
        pass
    return {}

# ══════════════════════════════════════════════════════════════════════════════
#  PHYSICAL MIGRATION MODEL
# ══════════════════════════════════════════════════════════════════════════════
def compute_band_position(bp: float, gel_pct: float, voltage: float,
                          buffer: str, time_h: float) -> float:
    """
    Return normalized gel position ∈ [0, 0.95].
      0   = at well (top)
      0.95 = near bottom of gel

    Physical model (Ferguson-Ogston inspired):
    ─────────────────────────────────────────
    • Sieving exponent α increases with gel%  → better resolution of small fragments
    • Speed ∝ voltage^0.92 × buffer_speed / gel%^0.40
    • Calibration: 1 kb at 100 V, 1% TAE, 1 h  →  position ≈ 0.40
    """
    if bp <= 0 or time_h <= 0:
        return 0.0

    bp       = float(max(bp, 1))
    gel_pct  = float(max(gel_pct, 0.1))
    voltage  = float(max(voltage, 1))

    alpha      = 0.55 + 0.125 * gel_pct          # e.g. 0.675 @ 1%, 0.80 @ 2%
    k_calib    = 0.40 * (1000.0 ** alpha)         # keep 1 kb standard at 0.40
    buf_speed  = BUF_FACTORS.get(buffer, BUF_FACTORS["TAE"])["speed"]
    gel_speed  = 1.0 / (gel_pct ** 0.40)
    v_factor   = (voltage / 100.0) ** 0.92

    pos = k_calib * time_h * v_factor * buf_speed * gel_speed / (bp ** alpha)
    return float(min(0.95, max(0.0, pos)))


def pos_to_y(pos: float) -> float:
    """Normalize position → matplotlib y-coordinate (0.97 = well, 0.04 = bottom)."""
    return 0.97 - 0.93 * pos

# ══════════════════════════════════════════════════════════════════════════════
#  PLASMID TOPOLOGY FORMS  (for plasmid isolation samples)
# ══════════════════════════════════════════════════════════════════════════════
# When DNA is extracted by miniprep, three topological forms co-migrate:
#
#  Form I  – Supercoiled (CCC): compact, migrates FASTEST.
#            Apparent size ≈ 0.60–0.70× actual bp (gel-% dependent).
#  Form II – Nicked / Relaxed circular (OC): migrates SLOWEST.
#            Apparent size ≈ 1.80–2.10× actual bp.
#  Form III– Linear: runs at true bp size (reference band).
#
# These ratios are well-established; supercoiling factor increases slightly
# with higher gel %, because compact forms are more differentially sieved.

def _topo_factors(gel_pct: float) -> tuple:
    """Return (sc_factor, oc_factor) apparent-size multipliers for the given gel %."""
    # Supercoiled runs faster → apparent bp is SMALLER than actual
    sc = max(0.58, 0.66 - 0.04 * (gel_pct - 1.0))   # 0.66 @ 1%, 0.58 @ 3%
    # Open-circle runs slower → apparent bp is LARGER than actual
    oc = min(2.20, 1.85 + 0.10 * (gel_pct - 1.0))   # 1.85 @ 1%, 2.05 @ 3%
    return sc, oc


def get_plasmid_forms(actual_bp: int, gel_pct: float) -> list:
    """
    Return list of (apparent_bp, label, alpha) tuples for the three forms.
    alpha = band opacity (supercoiled is brightest as most abundant in miniprep).
    """
    sc_f, oc_f = _topo_factors(gel_pct)
    return [
        (int(actual_bp * sc_f),  "SC",     1.00),   # Form I  – supercoiled
        (actual_bp,               "Lin",    0.65),   # Form III – linear
        (int(actual_bp * oc_f),  "OC",     0.75),   # Form II  – open circle
    ]

# ══════════════════════════════════════════════════════════════════════════════
#  RESTRICTION ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
def _find_enzymes(names):
    normed = [n.strip().replace("1","I").replace("2","II").replace("3","III")
              for n in names if n.strip()]
    return [e for e in AllEnzymes
            if any(getattr(e,"__name__","").lower() == n.lower() for n in normed)]


def compute_fragments(sequence: str, enzyme_names: list,
                      circular: bool = False, no_digest: bool = False) -> list:
    seq = _clean_seq(sequence)
    if not seq:
        return []
    if no_digest:
        return [len(seq)]

    s       = Seq(seq)
    enzymes = _find_enzymes(enzyme_names)
    if not enzymes:
        return [len(s)]

    cuts = sorted(set(p for ps in Analysis(enzymes, s).full().values() for p in ps))
    if not cuts:
        return [len(s)]

    if circular:
        if len(cuts) == 1:
            return [len(s)]
        sizes = []
        for i in range(len(cuts)):
            a, b = cuts[i], cuts[(i+1) % len(cuts)]
            sizes.append((len(s) - a + b) if b <= a else (b - a))
        return sizes

    spans = [0] + cuts + [len(s)]
    return [int(spans[i+1]-spans[i]) for i in range(len(spans)-1) if spans[i+1] > spans[i]]


def ladder_from_selection(selection: str, custom_str: str = "") -> list:
    if selection in PRESET_LADDERS and PRESET_LADDERS[selection]:
        return PRESET_LADDERS[selection]
    try:
        return sorted([int(float(p)) for p in custom_str.split(",") if p.strip()], reverse=True)
    except Exception:
        return []

# ══════════════════════════════════════════════════════════════════════════════
#  GEL RENDERING
# ══════════════════════════════════════════════════════════════════════════════
_BAND_COLOR   = (0.93, 0.97, 0.88)   # warm white / gel band
_LADDER_COLOR = (0.96, 0.96, 0.86)


def render_gel(ladder, samples_frags, samples_labels,
               gel_pct=1.0, voltage=100, buffer="TAE", time_h=1.0,
               dpi=150, img_w=700, img_h=1400, out_path=None):
    """
    Render an agarose gel image in SnapGene style:
      • Light grey outer background
      • Dark charcoal gel slab
      • Crisp thin white/grey bands
      • Size labels + right-pointing ticks on the left outside the gel
      • Lane numbers above the gel, no sample names inside
    Returns the matplotlib Figure.
    """
    # ── Figure / axes setup ───────────────────────────────────────────────────
    # Outer background = light grey (like SnapGene paper)
    _BG   = (0.91, 0.91, 0.91)   # light grey background
    _GEL  = (0.13, 0.13, 0.16)   # dark gel slab
    _SBAND = (1.00, 1.00, 1.00)  # sample bands: white
    _LBAND = (0.68, 0.68, 0.62)  # ladder bands: medium warm grey
    _LABEL = (0.10, 0.10, 0.10)  # size label text: near-black
    _TICK  = (0.30, 0.30, 0.30)  # tick colour

    fig = plt.figure(figsize=(img_w/dpi, img_h/dpi), dpi=dpi,
                     facecolor=_BG)

    # Axes occupy the right 68% of the figure; left 32% is label space
    ax = fig.add_axes([0.32, 0.06, 0.66, 0.88])
    ax.set_facecolor(_GEL)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.axis("off")

    n       = min(len(samples_frags), 9)
    x_ldr   = 0.10                      # ladder lane centre (in axes coords)
    _gap    = 0.088                     # fixed lane pitch
    _start  = 0.24
    xs = [_start + i * _gap for i in range(max(n, 1))]
    if xs and xs[-1] > 0.97:
        scale = (0.97 - _start) / (xs[-1] - _start)
        xs = [_start + (x - _start) * scale for x in xs]

    lane_w  = 0.055     # visual band width

    # ── Gel slab border ───────────────────────────────────────────────────────
    ax.add_patch(patches.Rectangle(
        (0.0, 0.0), 1.0, 1.0,
        linewidth=0.8, edgecolor=(0.40, 0.40, 0.40),
        facecolor=_GEL, zorder=0))

    # ── Buffer / well zone at top ─────────────────────────────────────────────
    _well_y   = 0.970
    _well_h   = 0.014
    _well_w   = lane_w

    # Draw a thin dark well slot for each lane
    for xc in [x_ldr] + xs:
        ax.add_patch(patches.Rectangle(
            (xc - _well_w/2, _well_y - _well_h/2), _well_w, _well_h,
            facecolor=(0.06, 0.06, 0.08), linewidth=0.4,
            edgecolor=(0.50, 0.50, 0.50), zorder=5))

    # ── Lane numbers ABOVE gel ────────────────────────────────────────────────
    # "MW" above ladder, "1","2"… above sample lanes — drawn in figure coords
    def _lane_fig_x(ax_x):
        """Convert axes-data x to figure fraction."""
        return ax.get_position().x0 + ax_x * ax.get_position().width

    fig.text(_lane_fig_x(x_ldr), ax.get_position().y1 + 0.005,
             "MW", ha="center", va="bottom",
             fontsize=6.5, color=(0.15, 0.15, 0.15), fontweight="bold",
             fontfamily="DejaVu Sans")

    for i, sx in enumerate(xs):
        fig.text(_lane_fig_x(sx), ax.get_position().y1 + 0.005,
                 str(i + 1), ha="center", va="bottom",
                 fontsize=6.5, color=(0.15, 0.15, 0.15),
                 fontfamily="DejaVu Sans")

    # "bp" label top-left of label area
    fig.text(0.01, ax.get_position().y1 + 0.005,
             "bp", ha="left", va="bottom",
             fontsize=6.5, color=_LABEL, fontweight="bold",
             fontfamily="DejaVu Sans")

    # ── Ladder bands + outside size labels with ticks ─────────────────────────
    ldr_left  = x_ldr - lane_w / 2   # left edge of ladder band
    ldr_right = x_ldr + lane_w / 2   # right edge

    # We'll skip duplicate y positions that are too close together
    prev_y = None
    _min_y_gap = 0.018   # minimum normalised gap before skipping a label

    for size in sorted(ladder, reverse=True):
        pos = compute_band_position(size, gel_pct, voltage, buffer, time_h)
        y   = pos_to_y(pos)

        # Ladder band — slightly thinner line style like SnapGene
        bh = 0.006
        ax.add_patch(patches.Rectangle(
            (ldr_left, y - bh/2), lane_w, bh,
            facecolor=_LBAND, linewidth=0, zorder=8))

        # Only label if not too close to previous label (avoids crowding)
        if prev_y is not None and abs(y - prev_y) < _min_y_gap:
            prev_y = y
            continue
        prev_y = y

        # Size label: drawn in figure coordinates, right-aligned with a tick
        # ax.ylim = (0, 1.05) so normalise y by ymax before converting to figure fraction
        _ymax = 1.05
        fig_y = ax.get_position().y0 + (y / _ymax) * ax.get_position().height

        # Format size
        if size >= 1000:
            lbl = f"{size:,}" if size < 10000 else f"{size:,}"
        else:
            lbl = str(size)

        # Label text (right-aligned against the gel left edge)
        fig.text(0.295, fig_y, lbl,
                 ha="right", va="center",
                 fontsize=5.8, color=_LABEL,
                 fontfamily="DejaVu Sans")

        # Short tick line from label → gel left edge (drawn in figure coords)
        tick_line = plt.Line2D(
            [0.298, ax.get_position().x0],
            [fig_y, fig_y],
            transform=fig.transFigure,
            color=_TICK, linewidth=0.5, zorder=10)
        fig.add_artist(tick_line)

    # ── Sample bands ──────────────────────────────────────────────────────────
    for i, frags in enumerate(samples_frags[:n]):
        sx = xs[i]
        for item in frags:
            if isinstance(item, (list, tuple)) and len(item) == 3:
                size, _lbl_ignored, band_alpha = item
            else:
                size, band_alpha = item, 1.0

            pos = compute_band_position(size, gel_pct, voltage, buffer, time_h)
            y   = pos_to_y(pos)
            bh  = 0.007 if size >= 2000 else 0.005

            ax.add_patch(patches.Rectangle(
                (sx - lane_w/2, y - bh/2), lane_w, bh,
                facecolor=(*_SBAND, band_alpha), linewidth=0, zorder=7))

    # ── Output ────────────────────────────────────────────────────────────────
    if out_path:
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        fig.savefig(out_path, dpi=600, bbox_inches="tight",
                    facecolor=fig.get_facecolor())

    return fig

# ══════════════════════════════════════════════════════════════════════════════
#  UI HELPERS
# ══════════════════════════════════════════════════════════════════════════════
def _btn(text, color=_BTN, font_size="12sp", **kw):
    b = Button(text=text, font_size=font_size, **kw)
    b.background_normal  = ""
    b.background_color   = color
    b.color              = _TEXT
    return b


def _lbl(text, size_hint=(1,None), height=dp(20), color=_TEXT, font_size="12sp", **kw):
    return Label(text=text, size_hint=size_hint, height=height,
                 color=color, font_size=font_size, **kw)


def _inp(**kw):
    return TextInput(
        foreground_color=(0.94,0.94,0.96,1),
        background_color=(0.20,0.20,0.27,1),
        cursor_color=(0.35,0.70,0.98,1),
        **kw)


def _panel_bg(widget, radius=10):
    """Attach a rounded-rectangle canvas background to widget."""
    with widget.canvas.before:
        Color(*_PANEL)
        rect = RoundedRectangle(pos=widget.pos, size=widget.size, radius=[radius])
    widget.bind(pos=lambda *_: setattr(rect, "pos", widget.pos),
                size=lambda *_: setattr(rect, "size", widget.size))

# ══════════════════════════════════════════════════════════════════════════════
#  NATIVE OS FILE / FOLDER DIALOGS  (tkinter → Windows Explorer / Finder / GTK)
# ══════════════════════════════════════════════════════════════════════════════
_SEQ_FILETYPES = [
    ("Sequence files", "*.fa *.fasta *.fas *.gb *.gbk *.genbank *.dna *.txt"),
    ("FASTA",          "*.fa *.fasta *.fas"),
    ("GenBank",        "*.gb *.gbk *.genbank"),
    ("SnapGene DNA",   "*.dna"),
    ("Plain text",     "*.txt"),
    ("All files",      "*.*"),
]


def open_native_file_dialog(callback):
    """
    Open the native OS file picker in a background thread (so Kivy stays responsive).
    When the user confirms, schedule `callback(records_dict, filepath)` on Kivy's
    main thread via Clock.schedule_once.
    """
    def _run():
        root = tk.Tk()
        root.withdraw()                    # hide the empty Tk window
        root.attributes("-topmost", True)  # dialog appears in front of Kivy
        path = _filedialog.askopenfilename(
            parent=root,
            title="Select Sequence File",
            initialdir=ROOT,
            filetypes=_SEQ_FILETYPES,
        )
        root.destroy()
        if path:
            records = parse_sequence_file(path)
            Clock.schedule_once(lambda _dt: callback(records, path))

    threading.Thread(target=_run, daemon=True).start()


def open_native_folder_dialog(callback, start_path=None):
    """
    Open the native OS folder picker in a background thread.
    When confirmed, schedule `callback(folder_path)` on Kivy's main thread.
    """
    def _run():
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        folder = _filedialog.askdirectory(
            parent=root,
            title="Choose Export Folder",
            initialdir=start_path or ROOT,
            mustexist=True,
        )
        root.destroy()
        if folder:
            Clock.schedule_once(lambda _dt: callback(folder))

    threading.Thread(target=_run, daemon=True).start()

# ══════════════════════════════════════════════════════════════════════════════
#  SAMPLE LANE WIDGET
# ══════════════════════════════════════════════════════════════════════════════
class SampleBox(BoxLayout):
    def __init__(self, idx, remove_cb, status_cb, **kw):
        super().__init__(orientation="vertical", spacing=4,
                         padding=(8,6,8,6), size_hint_y=None,
                         height=dp(262), **kw)          # +27 dp for plasmid row
        self.idx       = idx
        self._remove   = remove_cb
        self._status   = status_cb
        _panel_bg(self, radius=8)

        # ── Header ────────────────────────────────────────────────────────
        hdr = BoxLayout(size_hint_y=None, height=dp(30), spacing=6)
        self.title = Label(text=f"[b]Lane {idx}[/b]", markup=True,
                           color=_ACCENT, font_size="13sp", halign="left",
                           size_hint_x=1)
        hdr.add_widget(self.title)

        imp = _btn("📂 Import", color=_BTN, size_hint_x=None, width=dp(88), font_size="11sp")
        imp.bind(on_press=lambda *_: open_native_file_dialog(self._on_imported))
        hdr.add_widget(imp)

        rm = _btn("✕", color=_DANGER, size_hint_x=None, width=dp(34), font_size="13sp")
        rm.bind(on_press=lambda *_: remove_cb(self))
        hdr.add_widget(rm)
        self.add_widget(hdr)

        # ── Sequence textarea ──────────────────────────────────────────────
        self.seq = _inp(multiline=True, size_hint_y=None, height=dp(70),
                        hint_text="Paste sequence or use 📂 Import …")
        self.add_widget(self.seq)

        # ── Enzyme row ────────────────────────────────────────────────────
        er = BoxLayout(size_hint_y=None, height=dp(30), spacing=6)
        er.add_widget(_lbl("Enzymes:", size_hint=(None,None),
                           height=dp(28), width=dp(68), font_size="11sp"))
        self.enzymes = _inp(multiline=False, size_hint_x=1,
                            hint_text="EcoRI, BamHI …")
        er.add_widget(self.enzymes)
        self.pick = Spinner(text="＋ Add", values=COMMON_ENZYMES,
                            size_hint_x=None, width=dp(110), font_size="11sp")
        self.pick.bind(text=self._add_enz)
        er.add_widget(self.pick)
        self.add_widget(er)

        # ── Options row ───────────────────────────────────────────────────
        opt = BoxLayout(size_hint_y=None, height=dp(28), spacing=4)
        self.circ = CheckBox(size_hint_x=None, width=dp(22), color=_ACCENT)
        opt.add_widget(self.circ)
        opt.add_widget(_lbl("Circular", size_hint=(None,None),
                            height=dp(26), width=dp(60), font_size="11sp"))
        self.nodig = CheckBox(size_hint_x=None, width=dp(22), color=_ACCENT)
        opt.add_widget(self.nodig)
        opt.add_widget(_lbl("No digest", size_hint=(None,None),
                            height=dp(26), width=dp(72), font_size="11sp"))
        opt.add_widget(_lbl("Label:", size_hint=(None,None),
                            height=dp(26), width=dp(48), font_size="11sp"))
        self.lane_lbl = _inp(multiline=False, size_hint_x=1,
                             hint_text="optional")
        opt.add_widget(self.lane_lbl)
        self.add_widget(opt)

        # ── Plasmid isolation row ──────────────────────────────────────────
        plrow = BoxLayout(size_hint_y=None, height=dp(26), spacing=4)
        self.plasmid_cb = CheckBox(size_hint_x=None, width=dp(22), color=(0.98, 0.70, 0.22, 1))
        plrow.add_widget(self.plasmid_cb)
        plrow.add_widget(_lbl(
            "Plasmid isolation (shows SC / OC / Linear forms)",
            size_hint=(1, None), height=dp(24),
            color=(0.88, 0.72, 0.30, 1), font_size="10sp",
            halign="left", text_size=(None, None)))
        self.add_widget(plrow)

        # ── Fragment report ────────────────────────────────────────────────
        self.frag_lbl = _lbl("No fragments yet.", size_hint=(1,None),
                             height=dp(20), color=_HINT, font_size="10sp",
                             halign="left", text_size=(None,None))
        self.add_widget(self.frag_lbl)

    # ── public interface ──────────────────────────────────────────────────
    def set_index(self, i):
        self.idx = i
        self.title.text = f"[b]Lane {i}[/b]"

    def get_fragments(self, gel_pct: float = 1.0) -> list:
        seq  = _clean_seq(self.seq.text)
        if not seq:
            self.frag_lbl.text  = "No sequence entered."
            self.frag_lbl.color = _HINT
            return []
        enz   = [e.strip() for e in self.enzymes.text.split(",") if e.strip()]
        frags = sorted(compute_fragments(seq, enz,
                                        circular=self.circ.active,
                                        no_digest=self.nodig.active), reverse=True)

        # ── Plasmid isolation mode ─────────────────────────────────────────
        # If checked, replace each fragment with its three topological forms.
        # Only applies when no enzymes cut (undigested plasmid) — after
        # restriction digestion the fragments are linear anyway.
        enzymes_typed = bool([e for e in self.enzymes.text.split(",") if e.strip()])
        if self.plasmid_cb.active and not enzymes_typed and not self.nodig.active:
            topo_bands = []
            for bp in frags:
                for form in get_plasmid_forms(bp, gel_pct):
                    topo_bands.append(form)   # (apparent_bp, label, alpha)
            if topo_bands:
                form_sizes = [b[0] for b in topo_bands]
                lbl_parts  = [f"{b[1]}:{b[0]//1000:.1f}kb" if b[0]>=1000
                               else f"{b[1]}:{b[0]}bp" for b in topo_bands]
                self.frag_lbl.text  = "Forms: " + ", ".join(lbl_parts[:6])
                self.frag_lbl.color = (0.98, 0.82, 0.30, 1)
            return topo_bands

        if frags:
            parts = [f"{f//1000:.1f}kb" if f>=1000 else f"{f}bp" for f in frags[:7]]
            trail = f" … +{len(frags)-7}" if len(frags)>7 else ""
            self.frag_lbl.text  = "Fragments: " + ", ".join(parts) + trail
            self.frag_lbl.color = _GREEN
        else:
            self.frag_lbl.text  = "No cut sites found."
            self.frag_lbl.color = _HINT
        return frags

    # ── private ───────────────────────────────────────────────────────────
    def _add_enz(self, spinner, text):
        if text and text != "＋ Add":
            cur = self.enzymes.text.strip()
            self.enzymes.text = f"{cur}, {text}".lstrip(", ") if cur else text
            spinner.text = "＋ Add"

    def _on_imported(self, records, filepath):
        if not records:
            self._status(f"⚠ Could not parse: {os.path.basename(filepath)}", err=True)
            return
        name, seq = next(iter(records.items()))
        self.seq.text = seq
        if not self.lane_lbl.text:
            self.lane_lbl.text = name[:18]
        self._status(f"✓ Imported '{name}'  ({len(seq):,} bp)  ←  {os.path.basename(filepath)}")

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN APPLICATION
# ══════════════════════════════════════════════════════════════════════════════
class GelElectrophoresisApp(App):

    def build(self):
        self.title = "Gel Electrophoresis Simulator v2.0"
        self._sim_event    = None
        self._sim_t        = 0.1
        self._sim_run      = False
        self._debounce     = None
        self._preview_px   = os.path.join(tempfile.gettempdir(), "gel_preview.png")
        self._export_folder = os.path.expanduser("~")   # default export destination

        # ═══════════════════ ROOT ═════════════════════════════════════════
        root = BoxLayout(orientation="vertical", spacing=0, padding=0)

        # ── Top bar ───────────────────────────────────────────────────────
        tb = BoxLayout(size_hint_y=None, height=dp(46), padding=(12, 8))
        with tb.canvas.before:
            Color(0.09, 0.09, 0.12, 1); self._tb_bg = Rectangle()
        tb.bind(pos=lambda *_: setattr(self._tb_bg, "pos", tb.pos),
                size=lambda *_: setattr(self._tb_bg, "size", tb.size))
        tb.add_widget(Label(
            text="[b]🧬  Gel Electrophoresis Simulator  v2.0[/b]",
            markup=True, color=_ACCENT, font_size="16sp", halign="left"))
        root.add_widget(tb)

        # ══════════════════════════════════════════════════════════════════
        # BODY  →  3 columns:  [Settings 0.22] [Gel Image 0.44] [Lanes 0.34]
        # ══════════════════════════════════════════════════════════════════
        body = BoxLayout(orientation="horizontal", spacing=8, padding=(8, 6, 8, 4))

        # ── LEFT: Settings panel ──────────────────────────────────────────
        left = BoxLayout(orientation="vertical", spacing=6, size_hint_x=0.22)

        sc = BoxLayout(orientation="vertical", spacing=5, padding=(10, 8),
                       size_hint=(1, 1))
        _panel_bg(sc)

        sc.add_widget(Label(text="[b]⚙  Gel Settings[/b]", markup=True,
                            color=_ACCENT, font_size="12sp",
                            size_hint_y=None, height=dp(22), halign="left"))

        # Ladder spinner
        sc.add_widget(_lbl("Ladder:", size_hint=(1, None), height=dp(18),
                           color=_HINT, font_size="10sp"))
        self.ldr_spin = Spinner(text=list(PRESET_LADDERS)[0],
                                values=list(PRESET_LADDERS),
                                font_size="10sp", size_hint_y=None, height=dp(30))
        sc.add_widget(self.ldr_spin)

        self.custom_ldr = _inp(multiline=False, size_hint_y=None, height=dp(26),
                               hint_text="Custom ladder (bp, comma-sep)",
                               font_size="10sp")
        sc.add_widget(self.custom_ldr)

        # Gel % + Buffer on separate rows (narrow column)
        sc.add_widget(_lbl("Gel %:", size_hint=(1, None), height=dp(18),
                           color=_HINT, font_size="10sp"))
        self.gel_spin = Spinner(text="1.0%", values=GEL_OPTS,
                                font_size="10sp", size_hint_y=None, height=dp(30))
        sc.add_widget(self.gel_spin)

        sc.add_widget(_lbl("Buffer:", size_hint=(1, None), height=dp(18),
                           color=_HINT, font_size="10sp"))
        self.buf_spin = Spinner(text="TAE", values=BUF_OPTS,
                                font_size="10sp", size_hint_y=None, height=dp(30))
        sc.add_widget(self.buf_spin)

        self.buf_desc = _lbl(self._buf_hint("TAE"),
                             size_hint=(1, None), height=dp(28),
                             color=_HINT, font_size="9sp", halign="left",
                             text_size=(None, None))
        self.buf_spin.bind(text=lambda _, t: setattr(self.buf_desc, "text",
                                                      self._buf_hint(t)))
        sc.add_widget(self.buf_desc)

        # Voltage — slider + manual text input (kept in sync)
        sc.add_widget(_lbl("Voltage (V):", size_hint=(1, None), height=dp(18),
                           color=_HINT, font_size="10sp"))
        vrow = BoxLayout(size_hint_y=None, height=dp(30), spacing=4)
        self.v_slider = Slider(min=30, max=300, value=100, step=1, size_hint_x=1)
        self.v_inp = _inp(multiline=False, size_hint_x=None, width=dp(52),
                          text="100", font_size="11sp",
                          input_filter="int", padding=(4, 6))
        # slider → text
        self.v_slider.bind(value=lambda _, v: setattr(self.v_inp, "text", str(int(v))))
        self.v_slider.bind(value=self._defer_preview)
        # text → slider (on Enter or focus-out)
        self.v_inp.bind(on_text_validate=self._sync_v_inp)
        self.v_inp.bind(focus=lambda w, f: self._sync_v_inp(w) if not f else None)
        vrow.add_widget(self.v_slider)
        vrow.add_widget(self.v_inp)
        sc.add_widget(vrow)

        # Run time — slider + manual text input (kept in sync)
        sc.add_widget(_lbl("Run time (h):", size_hint=(1, None), height=dp(18),
                           color=_HINT, font_size="10sp"))
        trow = BoxLayout(size_hint_y=None, height=dp(30), spacing=4)
        self.t_slider = Slider(min=0.05, max=6.0, value=1.0, step=0.05, size_hint_x=1)
        self.t_inp = _inp(multiline=False, size_hint_x=None, width=dp(52),
                          text="1.00", font_size="11sp", padding=(4, 6))
        # slider → text
        self.t_slider.bind(value=lambda _, v: setattr(self.t_inp, "text", f"{v:.2f}"))
        self.t_slider.bind(value=self._defer_preview)
        # text → slider
        self.t_inp.bind(on_text_validate=self._sync_t_inp)
        self.t_inp.bind(focus=lambda w, f: self._sync_t_inp(w) if not f else None)
        trow.add_widget(self.t_slider)
        trow.add_widget(self.t_inp)
        sc.add_widget(trow)

        # Spacer
        sc.add_widget(Label(size_hint_y=1))

        # Simulate / Reset buttons
        self.play_btn = _btn("▶  Simulate", color=_GREEN,
                             font_size="11sp", size_hint_y=None, height=dp(34))
        self.play_btn.bind(on_press=self._toggle_sim)
        sc.add_widget(self.play_btn)

        self.stop_btn = _btn("⏮  Reset", font_size="11sp",
                             size_hint_y=None, height=dp(30))
        self.stop_btn.bind(on_press=self._reset_sim)
        sc.add_widget(self.stop_btn)

        # Export folder chooser
        sc.add_widget(_lbl("Export folder:", size_hint=(1, None), height=dp(18),
                           color=_HINT, font_size="10sp"))
        self.export_lbl = _lbl(
            self._short_path(self._export_folder),
            size_hint=(1, None), height=dp(24),
            color=_TEXT, font_size="9sp", halign="left",
            text_size=(None, None))
        sc.add_widget(self.export_lbl)

        folder_btn = _btn("📁  Choose Folder", font_size="10sp",
                          size_hint_y=None, height=dp(30))
        folder_btn.bind(on_press=lambda *_: open_native_folder_dialog(
            self._on_folder_chosen,
            start_path=self._export_folder))
        sc.add_widget(folder_btn)

        exp_btn = _btn("💾  Export 600 DPI PNG", color=(0.28, 0.28, 0.50, 1),
                       font_size="11sp", size_hint_y=None, height=dp(34))
        exp_btn.bind(on_press=self._export)
        sc.add_widget(exp_btn)

        left.add_widget(sc)
        body.add_widget(left)

        # ── CENTER: Gel image (dominant) ──────────────────────────────────
        center = BoxLayout(orientation="vertical", spacing=4, size_hint_x=0.44)

        self.gel_img = Image(size_hint=(1, 1), allow_stretch=True, keep_ratio=True)
        _panel_bg(self.gel_img, radius=8)
        center.add_widget(self.gel_img)

        analyze_btn = _btn("🔬  Analyze & Preview Gel",
                           color=(0.25, 0.50, 0.28, 1),
                           font_size="12sp", size_hint_y=None, height=dp(38))
        analyze_btn.bind(on_press=lambda *_: self._analyze())
        center.add_widget(analyze_btn)

        body.add_widget(center)

        # ── RIGHT: Sample lanes ───────────────────────────────────────────
        right = BoxLayout(orientation="vertical", spacing=6, size_hint_x=0.34)

        rh = BoxLayout(size_hint_y=None, height=dp(40), spacing=8, padding=(0, 4))
        rh.add_widget(Label(text="[b]🧪  Sample Lanes[/b]", markup=True,
                            color=_ACCENT, font_size="13sp"))
        add_b = _btn("＋  Lane", color=_ACCENT,
                     size_hint_x=None, width=dp(88), font_size="11sp")
        add_b.bind(on_press=lambda *_: self._add_lane())
        rh.add_widget(add_b)
        right.add_widget(rh)

        self.lane_grid = GridLayout(cols=1, spacing=8, size_hint_y=None, padding=(0, 4))
        self.lane_grid.bind(minimum_height=self.lane_grid.setter("height"))
        self.lanes: list = []
        sv = ScrollView(size_hint=(1, 1))
        sv.add_widget(self.lane_grid)
        right.add_widget(sv)
        body.add_widget(right)

        root.add_widget(body)

        # ── Status bar ────────────────────────────────────────────────────
        sb = BoxLayout(size_hint_y=None, height=dp(28), padding=(10, 4))
        with sb.canvas.before:
            Color(0.07, 0.07, 0.09, 1); self._sb_bg = Rectangle()
        sb.bind(pos=lambda *_: setattr(self._sb_bg, "pos", sb.pos),
                size=lambda *_: setattr(self._sb_bg, "size", sb.size))
        self.status = Label(
            text="Ready — add sample lanes and click  🔬 Analyze & Preview.",
            color=_HINT, font_size="10sp", halign="left", text_size=(None, None))
        sb.add_widget(self.status)
        root.add_widget(sb)

        # Initial lanes + first render
        for _ in range(2):
            self._add_lane()
        Clock.schedule_once(lambda *_: self._analyze(), 0.6)
        return root

    # ── Helpers ───────────────────────────────────────────────────────────
    def _short_path(self, path, max_len=28):
        return path if len(path) <= max_len else "…" + path[-(max_len - 1):]

    def _buf_hint(self, buf):
        hints = {
            "TAE": "TAE — faster migration, ideal for large fragments (≥2 kb)",
            "TBE": "TBE — better resolution, ideal for small fragments (≤1 kb)",
        }
        return hints.get(buf, "")

    def _set_status(self, text, err=False):
        self.status.text  = text
        self.status.color = _DANGER if err else (0.45, 0.88, 0.58, 1)

    def _on_folder_chosen(self, folder):
        self._export_folder = folder
        self.export_lbl.text = self._short_path(folder)
        self._set_status(f"Export folder set →  {folder}")

    def _get_params(self):
        gel_pct = float(self.gel_spin.text.strip("%"))
        # Read voltage and time from text inputs (authoritative); slider is display-only
        try:
            voltage = max(30, min(300, int(self.v_inp.text.strip())))
        except ValueError:
            voltage = int(self.v_slider.value)
        try:
            time_h = max(0.05, min(6.0, round(float(self.t_inp.text.strip()), 3)))
        except ValueError:
            time_h = round(self.t_slider.value, 3)
        buffer  = self.buf_spin.text
        ladder  = ladder_from_selection(self.ldr_spin.text,
                                        self.custom_ldr.text.strip())
        return gel_pct, voltage, buffer, time_h, ladder

    # ── Slider ↔ text-input sync ──────────────────────────────────────────
    def _sync_v_inp(self, widget):
        """Text input → slider for voltage."""
        try:
            v = max(30, min(300, int(widget.text.strip())))
        except ValueError:
            v = int(self.v_slider.value)
        widget.text = str(v)
        self.v_slider.value = v
        self._defer_preview()

    def _sync_t_inp(self, widget):
        """Text input → slider for run time."""
        try:
            t = max(0.05, min(6.0, round(float(widget.text.strip()), 2)))
        except ValueError:
            t = round(self.t_slider.value, 2)
        widget.text = f"{t:.2f}"
        self.t_slider.value = t
        self._defer_preview()

    # ── Lane management ───────────────────────────────────────────────────
    def _add_lane(self):
        if len(self.lanes) >= 9:
            self._set_status("⚠  Maximum 9 sample lanes reached.", err=True)
            return
        b = SampleBox(len(self.lanes)+1, self._remove_lane, self._set_status)
        self.lanes.append(b)
        self.lane_grid.add_widget(b)

    def _remove_lane(self, box):
        if box in self.lanes:
            self.lane_grid.remove_widget(box)
            self.lanes.remove(box)
            for i, b in enumerate(self.lanes, 1):
                b.set_index(i)

    # ── Analysis / render ─────────────────────────────────────────────────
    def _analyze(self):
        gel_pct, voltage, buffer, time_h, ladder = self._get_params()
        if not ladder:
            self._set_status("⚠  Select a ladder or enter custom sizes.", err=True)
            return

        frags_list = [sb.get_fragments(gel_pct=gel_pct) for sb in self.lanes]
        labels     = [sb.lane_lbl.text.strip() for sb in self.lanes]

        try:
            fig = render_gel(ladder, frags_list, labels,
                             gel_pct=gel_pct, voltage=voltage,
                             buffer=buffer, time_h=time_h,
                             dpi=110, img_w=520, img_h=1100)
            fig.savefig(self._preview_px, dpi=110, bbox_inches="tight",
                        facecolor=fig.get_facecolor())
            plt.close(fig)
            self.gel_img.source = self._preview_px
            self.gel_img.reload()
            self._set_status(
                f"Preview  ·  {gel_pct}% gel  |  {buffer}  |  {int(voltage)} V  |  t={time_h:.2f} h")
        except Exception as exc:
            self._set_status(f"⚠  Render error: {exc}", err=True)

    def _defer_preview(self, *_):
        if self._debounce:
            self._debounce.cancel()
        self._debounce = Clock.schedule_once(lambda *_: self._analyze(), 0.30)

    # ── Simulation (animated time slider) ────────────────────────────────
    def _toggle_sim(self, *_):
        if self._sim_run:
            self._pause_sim()
        else:
            self._start_sim()

    def _start_sim(self):
        # Determine target run time from user input (t_inp), clamped to slider range
        try:
            target_h = max(self.t_slider.min,
                           min(self.t_slider.max,
                               round(float(self.t_inp.text.strip()), 2)))
        except ValueError:
            target_h = self.t_slider.max

        self._sim_run    = True
        self._sim_target = target_h
        # Always animate from the beginning (t = slider.min)
        self._sim_t      = self.t_slider.min
        self.t_slider.value = self._sim_t
        self.t_inp.text     = f"{self._sim_t:.2f}"

        self.play_btn.text = "⏸  Pause"
        self.play_btn.background_color = (0.75, 0.55, 0.15, 1)
        self._set_status(f"▶ Running protocol: {int(self.v_slider.value)} V "
                         f"for {target_h:.2f} h …")

        def step(dt):
            self._sim_t = round(self._sim_t + 0.05, 3)
            if self._sim_t >= self._sim_target:
                self._sim_t = self._sim_target
                self.t_slider.value = self._sim_t
                self.t_inp.text     = f"{self._sim_t:.2f}"
                self._stop_sim()
                self._set_status(
                    f"✓ Run complete — {int(self.v_slider.value)} V  ×  {self._sim_target:.2f} h")
                return False
            self.t_slider.value = self._sim_t
            self.t_inp.text     = f"{self._sim_t:.2f}"

        self._sim_event = Clock.schedule_interval(step, 0.14)

    def _pause_sim(self):
        self._sim_run = False
        self.play_btn.text = "▶  Simulate Run"
        self.play_btn.background_color = _GREEN
        if self._sim_event:
            self._sim_event.cancel()
            self._sim_event = None

    def _stop_sim(self):
        self._pause_sim()

    def _reset_sim(self, *_):
        self._stop_sim()
        self.t_slider.value = self.t_slider.min
        self._analyze()

    # ── Export ───────────────────────────────────────────────────────────
    def _export(self, *_):
        gel_pct, voltage, buffer, time_h, ladder = self._get_params()
        if not ladder:
            self._set_status("⚠  No ladder selected.", err=True)
            return

        frags_list = [sb.get_fragments(gel_pct=gel_pct) for sb in self.lanes]
        labels     = [sb.lane_lbl.text.strip() for sb in self.lanes]
        stamp      = datetime.now().strftime("%Y%m%d_%H%M%S")
        out        = os.path.join(self._export_folder, f"gel_{stamp}.png")

        try:
            fig = render_gel(ladder, frags_list, labels,
                             gel_pct=gel_pct, voltage=voltage,
                             buffer=buffer, time_h=time_h,
                             dpi=600, img_w=919, img_h=2214,
                             out_path=out)
            plt.close(fig)
            self._set_status(f"✓ Exported →  {out}")
        except Exception as exc:
            self._set_status(f"⚠  Export error: {exc}", err=True)


# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    GelElectrophoresisApp().run()
