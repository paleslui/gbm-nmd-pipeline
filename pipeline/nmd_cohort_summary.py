#!/usr/bin/env python3
"""
nmd_cohort_summary.py — GBM Pipeline Stage 3: NMD Cohort Summary
Part of the GBM NMD-Neoantigen Pipeline.
https://github.com/paleslui/gbm-nmd-pipeline

Aggregates per-sample NMD scoring outputs into a single cohort report.
Focuses on NMD-actionable neoantigen candidates (frameshift only).
Missense, inframe and stop_gained variants are summarized at the top
(variant landscape) but excluded from downstream NMD analysis.

Outputs (in --out_dir):
  cohort_candidates.tsv     all candidates, all variant types
                            (now also: binder_class, clonality columns)
  cohort_summary.tsv        cohort-level counts (variant landscape + FS NMD)
  cohort_paired.tsv         per-patient T (primary) vs M (recurrent), FS only
  cohort_tier1.tsv          TIER1 high-priority candidates only
  cohort_landscape.tsv      variant-type landscape

  per_sample_nmd_summary.tsv             Module 3: variant-level NMD aggregation
  per_sample_neoepitope_summary.tsv      Module 4: binder/clonality burden
  per_sample_neoepitope_nmd_summary.tsv  Module 5: neoepitopes split by NMD class
  paired_stats_stage3.tsv                Wilcoxon + BH-FDR for Modules 3/4/5
  m3_*.png m4_*.png m5_*.png             12 module plots (300 dpi)
  cohort_report.html        single HTML report (intro + cohort overview +
                            landscape + Modules 3/4/5 + per-pair overlap widget +
                            per-patient drill-down widget + genes/HLA/tiers)

Optional Stage 1 join inputs:
  --stage1_summary <sample_mutation_summary.tsv>  for truncating_total / TMB
  --paired_overlap <paired_variant_overlap.tsv>   for the overlap widget
  Both degrade gracefully (NaN-filled / skipped) when absent.

Usage:
  python nmd_cohort_summary.py --input_dir <per_sample_dir> --out_dir <cohort_dir> \\
      [--stage1_summary <tsv>] [--paired_overlap <tsv>]
"""
import argparse, base64, io, json, re, sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, wilcoxon

# Reuse plotting + style from per-sample scorer for consistent look
sys.path.insert(0, str(Path(__file__).parent))
from nmd_scoring import (
    plot_tiers, plot_ic50, plot_nmd_breakdown, plot_confidence,
    hla_allele_breakdown, _b64fig, COL_S, COL_I, COL_U,
)

SAMPLE_RE = re.compile(r"^(\d+)_([TM])$")

# Variant types kept for NMD-actionable downstream analysis
NMD_ACTIONABLE = {"FS"}

# ─── Cohort-specific palette (extra colors not in nmd_scoring.py) ────────────
COL_T = "#378ADD"   # primary
COL_M = "#D85A30"   # recurrent
COL_TIER1 = "#F0C040"
COL_TIER2 = "#F5944D"
COL_TIER3 = "#378ADD"
COL_FS = "#D85A30"
COL_MIS = "#888888"
COL_INFRAME = "#bbbbbb"


# ═════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═════════════════════════════════════════════════════════════════════════════

def load_cohort(per_sample_dir: Path) -> pd.DataFrame:
    """Load per-sample TSVs from per_sample_dir/<sample>/nmd_scored_candidates.tsv"""
    rows = []
    sample_dirs = sorted(p for p in per_sample_dir.iterdir() if p.is_dir())
    print(f"[INFO] Found {len(sample_dirs)} per-sample dirs")
    n_loaded, n_skipped = 0, 0
    for d in sample_dirs:
        m = SAMPLE_RE.match(d.name)
        if not m:
            n_skipped += 1
            continue
        patient, tp = m.group(1), m.group(2)
        sample = f"{patient}_{tp}"
        timepoint = "primary" if tp == "T" else "recurrent"
        tsv = d / "nmd_scored_candidates.tsv"
        if not tsv.is_file():
            n_skipped += 1
            continue
        try:
            df = pd.read_csv(tsv, sep="\t", dtype=str)
        except pd.errors.EmptyDataError:
            n_skipped += 1
            continue
        if df.empty:
            n_skipped += 1
            continue
        df.insert(0, "sample", sample)
        df.insert(1, "patient", patient)
        df.insert(2, "timepoint", timepoint)
        rows.append(df)
        n_loaded += 1
    print(f"[INFO] Loaded {n_loaded} samples, skipped {n_skipped}")
    if not rows:
        return pd.DataFrame()
    cohort = pd.concat(rows, ignore_index=True)
    for col in ("best_ic50", "median_ic50", "nmd_confidence_score"):
        if col in cohort.columns:
            cohort[col] = pd.to_numeric(cohort[col], errors="coerce")
    print(f"[INFO] Cohort: {len(cohort)} candidates from {cohort['sample'].nunique()} samples")
    return cohort


# ═════════════════════════════════════════════════════════════════════════════
# AGGREGATE TABLES
# ═════════════════════════════════════════════════════════════════════════════

def variant_landscape(cohort: pd.DataFrame) -> pd.DataFrame:
    """Variant type breakdown across all candidates (the 'landscape')."""
    if cohort.empty:
        return pd.DataFrame()
    counts = cohort["Variant Type"].value_counts()
    rows = [{
        "Variant Type": vt,
        "Count": int(counts.get(vt, 0)),
        "% of total": f"{100*counts.get(vt, 0)/len(cohort):.1f}",
        "NMD-actionable?": "yes" if vt in NMD_ACTIONABLE else "no"
    } for vt in counts.index]
    return pd.DataFrame(rows)


def cohort_summary(cohort: pd.DataFrame, fs_cohort: pd.DataFrame, n_input_samples: int) -> pd.DataFrame:
    """Cohort-level aggregate counts for both the full landscape and the FS subset."""
    if cohort.empty:
        return pd.DataFrame([("No candidates", 0)], columns=["Metric", "Count"])
    rows = [
        ("─── Input ───", ""),
        ("Samples loaded",                        cohort["sample"].nunique()),
        ("Samples with 0 candidates",             n_input_samples - cohort["sample"].nunique()),
        ("Patients represented",                  cohort["patient"].nunique()),
        ("Total candidates (all variant types)",  len(cohort)),
        ("", ""),
        ("─── NMD-actionable (frameshift only) ───", ""),
        ("FS candidates",                         len(fs_cohort)),
        ("FS samples with ≥1 candidate",          fs_cohort["sample"].nunique() if len(fs_cohort) else 0),
        ("FS patients represented",               fs_cohort["patient"].nunique() if len(fs_cohort) else 0),
        ("", ""),
        ("NMD-SENSITIVE",                         int((fs_cohort["nmd_consensus"]=="SENSITIVE").sum())),
        ("NMD-INSENSITIVE",                       int((fs_cohort["nmd_consensus"]=="INSENSITIVE").sum())),
        ("UNCERTAIN / UNKNOWN",                   int(fs_cohort["nmd_consensus"].isin(["UNCERTAIN","UNKNOWN"]).sum())),
        ("", ""),
        ("TIER1 (NMD-sensitive, IC50<50)",        int((fs_cohort["priority_tier"]=="TIER1").sum())),
        ("TIER2 (NMD-sensitive, IC50<500)",       int((fs_cohort["priority_tier"]=="TIER2").sum())),
        ("TIER3 controls (NMD-insensitive)",      int((fs_cohort["priority_tier"]=="TIER3_control").sum())),
        ("Unclassified",                          int((fs_cohort["priority_tier"]=="UNCLASSIFIED").sum())),
        ("", ""),
        ("High confidence (3 — both methods agree)", int((fs_cohort["nmd_confidence_score"]==3).sum())),
        ("Medium confidence (2 — single method)",    int((fs_cohort["nmd_confidence_score"]==2).sum())),
        ("Low confidence (1 — methods disagree)",    int((fs_cohort["nmd_confidence_score"]==1).sum())),
        ("No data (0)",                              int((fs_cohort["nmd_confidence_score"]==0).sum())),
    ]
    return pd.DataFrame(rows, columns=["Metric", "Count"])


def per_patient_paired(fs_cohort: pd.DataFrame) -> pd.DataFrame:
    """T vs M comparison per patient on FS-only cohort."""
    if fs_cohort.empty:
        return pd.DataFrame()
    rows = []
    for patient, sub in fs_cohort.groupby("patient"):
        prim  = sub[sub["timepoint"] == "primary"]
        recur = sub[sub["timepoint"] == "recurrent"]
        rows.append({
            "patient":     patient,
            "T_total":     len(prim),
            "M_total":     len(recur),
            "T_tier1":     int((prim["priority_tier"] == "TIER1").sum()),
            "M_tier1":     int((recur["priority_tier"] == "TIER1").sum()),
            "T_sensitive": int((prim["nmd_consensus"] == "SENSITIVE").sum()),
            "M_sensitive": int((recur["nmd_consensus"] == "SENSITIVE").sum()),
            "delta_total": len(recur) - len(prim),
            "delta_tier1": int((recur["priority_tier"] == "TIER1").sum()
                               - (prim["priority_tier"] == "TIER1").sum()),
        })
    return pd.DataFrame(rows).sort_values(
        ["delta_tier1","delta_total"], ascending=[False,False]).reset_index(drop=True)


# ═════════════════════════════════════════════════════════════════════════════
# COHORT-SPECIFIC PLOTS
# ═════════════════════════════════════════════════════════════════════════════

def plot_variant_landscape(cohort: pd.DataFrame) -> str:
    """Bar chart of variant types across all candidates (the landscape)."""
    if cohort.empty:
        return ""
    counts = cohort["Variant Type"].value_counts()
    fig, ax = plt.subplots(figsize=(7, 3.8))
    colors = []
    for v in counts.index:
        if v == "FS": colors.append(COL_FS)
        elif v == "missense": colors.append(COL_MIS)
        else: colors.append(COL_INFRAME)
    bars = ax.bar(range(len(counts)), counts.values, color=colors)
    for i, (v, c) in enumerate(zip(counts.values, counts.index)):
        pct = 100 * v / len(cohort)
        ax.text(i, v + max(counts.values) * 0.02, f"{v}\n({pct:.1f}%)",
                ha="center", va="bottom", fontsize=10)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, fontsize=11)
    ax.set_ylabel("Candidates")
    ax.set_title(f"Variant landscape — {len(cohort)} candidates across cohort", fontsize=12)
    ax.set_ylim(top=max(counts.values) * 1.18)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    return _b64fig(fig)


def plot_tier_by_timepoint(fs_cohort: pd.DataFrame) -> str:
    if fs_cohort.empty:
        return ""
    counts = (fs_cohort.groupby(["timepoint","priority_tier"]).size()
              .unstack(fill_value=0)
              .reindex(columns=["TIER1","TIER2","TIER3_control","UNCLASSIFIED"], fill_value=0))
    fig, ax = plt.subplots(figsize=(9, 4))
    bar_w = 0.35
    x = range(len(counts.columns))
    primary   = counts.loc["primary"]   if "primary"   in counts.index else [0]*4
    recurrent = counts.loc["recurrent"] if "recurrent" in counts.index else [0]*4
    b1 = ax.bar([i - bar_w/2 for i in x], primary,   bar_w, color=COL_T, label="Primary (T)")
    b2 = ax.bar([i + bar_w/2 for i in x], recurrent, bar_w, color=COL_M, label="Recurrent (M)")
    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(bar.get_x() + bar.get_width()/2, h + 0.3, f"{int(h)}",
                        ha="center", fontsize=9)
    ax.set_xticks(list(x))
    ax.set_xticklabels(["TIER1","TIER2","TIER3 ctrl","Unclassified"], fontsize=10)
    ax.set_ylabel("Frameshift candidates")
    ax.set_title("Tier distribution: primary (T) vs recurrent (M) — FS only", fontsize=12)
    ax.legend(fontsize=9)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    return _b64fig(fig)


def plot_paired(paired: pd.DataFrame) -> str:
    if paired.empty:
        return ""
    fig, ax = plt.subplots(figsize=(11, 4.5))
    bar_w = 0.4
    x = range(len(paired))
    ax.bar([i - bar_w/2 for i in x], paired["T_tier1"], bar_w, color=COL_T, label="Primary (T)")
    ax.bar([i + bar_w/2 for i in x], paired["M_tier1"], bar_w, color=COL_M, label="Recurrent (M)")
    ax.set_xticks(list(x))
    ax.set_xticklabels(paired["patient"], rotation=90, fontsize=8)
    ax.set_xlabel("Patient")
    ax.set_ylabel("TIER1 candidates")
    ax.set_title("TIER1 NMD-sensitive neoantigens per patient: primary vs recurrent (FS only)",
                 fontsize=12)
    ax.legend(fontsize=9)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    return _b64fig(fig)


def plot_top_genes(fs_cohort: pd.DataFrame, n: int = 20) -> str:
    if fs_cohort.empty:
        return ""
    tier1 = fs_cohort[fs_cohort["priority_tier"] == "TIER1"]
    if tier1.empty:
        return ""
    counts = tier1["Gene Name"].value_counts().head(n)
    fig, ax = plt.subplots(figsize=(8, max(3, 0.32 * len(counts))))
    counts.iloc[::-1].plot(kind="barh", ax=ax, color=COL_TIER1)
    for i, v in enumerate(counts.iloc[::-1]):
        ax.text(v + 0.05, i, str(int(v)), va="center", fontsize=9)
    ax.set_xlabel("TIER1 candidate count")
    ax.set_title(f"Top {len(counts)} genes producing TIER1 candidates (FS only)", fontsize=12)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    return _b64fig(fig)


# ═════════════════════════════════════════════════════════════════════════════
# MAX MODULES 3/4/5 — variant-level NMD, neoepitope burden, neoepitopes-by-NMD
# ═════════════════════════════════════════════════════════════════════════════

# Variant identity = (chrom, pos, ref, alt, sample). NMD is a transcript-level
# property, so all peptides of one variant share a consensus; worst-case folding
# (SENSITIVE > INSENSITIVE > UNCERTAIN > UNKNOWN) is defensive.
_NMD_PRIORITY = ["SENSITIVE", "INSENSITIVE", "UNCERTAIN", "UNKNOWN"]
_NMD_RANK = {c: i for i, c in enumerate(_NMD_PRIORITY)}


def _variant_keys(df: pd.DataFrame) -> pd.Series:
    """Build a per-row variant identity key 'chrom:pos:ref>alt@sample'."""
    return (df["Chromosome"].astype(str) + ":" + df["Start"].astype(str) + ":"
            + df["Reference"].astype(str) + ">" + df["Variant"].astype(str)
            + "@" + df["sample"].astype(str))


def augment_candidates(cohort: pd.DataFrame) -> pd.DataFrame:
    """Append binder_class and clonality columns to the candidate table.

    Inputs:  cohort candidate DataFrame (must have 'Best MT IC50 Score',
             'Tumor DNA VAF').
    Output:  same DataFrame with two appended columns:
             binder_class — Strong (IC50<50 nM) / Weak (50–500) / Non-binder
                            (≥500 or null);
             clonality    — Clonal (VAF≥0.30) / Subclonal (<0.30) / Unknown (null).
    Side effects: mutates and returns ``cohort``.
    """
    if cohort.empty:
        cohort["binder_class"] = pd.Series(dtype=str)
        cohort["clonality"] = pd.Series(dtype=str)
        return cohort
    ic50 = pd.to_numeric(cohort["Best MT IC50 Score"], errors="coerce")
    vaf = pd.to_numeric(cohort["Tumor DNA VAF"], errors="coerce")

    def _binder(x: float) -> str:
        if pd.isna(x):
            return "Non-binder"
        if x < 50:
            return "Strong"
        if x < 500:
            return "Weak"
        return "Non-binder"

    def _clon(x: float) -> str:
        if pd.isna(x):
            return "Unknown"
        return "Clonal" if x >= 0.30 else "Subclonal"

    cohort["binder_class"] = ic50.map(_binder)
    cohort["clonality"] = vaf.map(_clon)
    return cohort


def variant_nmd_table(cohort: pd.DataFrame) -> pd.DataFrame:
    """Collapse the candidate table to one row per variant with worst-case NMD.

    Input:   cohort candidate DataFrame.
    Output:  DataFrame with columns vkey, sample, patient, timepoint,
             Variant Type, nmd (worst-case consensus across the variant's
             peptides). One row per unique (chrom,pos,ref,alt,sample).
    """
    if cohort.empty:
        return pd.DataFrame(columns=["vkey", "sample", "patient", "timepoint",
                                     "Variant Type", "nmd"])
    df = cohort[["sample", "patient", "timepoint", "Variant Type",
                 "nmd_consensus"]].copy()
    df["vkey"] = _variant_keys(cohort)
    df["nmd"] = df["nmd_consensus"].where(
        df["nmd_consensus"].isin(_NMD_PRIORITY), "UNKNOWN")
    df["_rank"] = df["nmd"].map(_NMD_RANK)
    # worst-case = highest priority = smallest rank
    df = df.sort_values("_rank").groupby("vkey", as_index=False).first()
    return df[["vkey", "sample", "patient", "timepoint", "Variant Type", "nmd"]]


def _sample_meta(cohort: pd.DataFrame) -> pd.DataFrame:
    """Unique (sample, patient, timepoint) rows, ordered patient↑ then T,M."""
    meta = cohort[["sample", "patient", "timepoint"]].drop_duplicates()
    meta = meta.sort_values(
        ["patient", "timepoint"],
        key=lambda s: s.map(int) if s.name == "patient"
        else s.map({"primary": 0, "recurrent": 1})).reset_index(drop=True)
    return meta


def module3_nmd_summary(cohort: pd.DataFrame, variant_tbl: pd.DataFrame,
                        stage1: Optional[pd.DataFrame]) -> pd.DataFrame:
    """Module 3 — per-sample variant-level NMD aggregation (FS / truncating).

    Inputs:  cohort table, the per-variant NMD table, and the optional Stage 1
             sample_mutation_summary (joined on 'sample' for truncating_total).
    Output:  one row per sample (all samples with ≥1 candidate) with columns
             sample, patient, timepoint, truncating_total_from_stage1,
             truncating_in_pvacseq, nmd_sensitive, nmd_insensitive,
             nmd_uncertain, nmd_unknown, nmd_fraction. Category counts are over
             unique FS variants; nmd_fraction = sensitive/(sensitive+insensitive),
             NaN if that denominator is zero.
    """
    meta = _sample_meta(cohort)
    fs = variant_tbl[variant_tbl["Variant Type"] == "FS"]
    s1 = (stage1.set_index("sample")["truncating_total"].to_dict()
          if stage1 is not None and "truncating_total" in stage1.columns else {})

    rows: List[Dict[str, object]] = []
    for _, m in meta.iterrows():
        sub = fs[fs["sample"] == m["sample"]]
        sens = int((sub["nmd"] == "SENSITIVE").sum())
        insens = int((sub["nmd"] == "INSENSITIVE").sum())
        unc = int((sub["nmd"] == "UNCERTAIN").sum())
        unk = int((sub["nmd"] == "UNKNOWN").sum())
        denom = sens + insens
        rows.append({
            "sample": m["sample"], "patient": m["patient"],
            "timepoint": m["timepoint"],
            "truncating_total_from_stage1": float(s1.get(m["sample"], np.nan)),
            "truncating_in_pvacseq": len(sub),
            "nmd_sensitive": sens, "nmd_insensitive": insens,
            "nmd_uncertain": unc, "nmd_unknown": unk,
            "nmd_fraction": round(sens / denom, 4) if denom else np.nan,
        })
    return pd.DataFrame(rows)


def module4_neoepitope_summary(cohort: pd.DataFrame) -> pd.DataFrame:
    """Module 4 — per-sample neoepitope (binder) burden and clonality.

    Input:   cohort table augmented with binder_class / clonality.
    Output:  one row per sample with columns sample, patient, timepoint,
             total_peptides, strong_binders, weak_binders, non_binders,
             clonal_binders, subclonal_binders, fraction_clonal. Counts are over
             all peptide rows; clonal/subclonal counts exclude Non-binders;
             fraction_clonal = clonal/(clonal+subclonal), NaN if denominator zero.
    """
    meta = _sample_meta(cohort)
    rows: List[Dict[str, object]] = []
    for _, m in meta.iterrows():
        sub = cohort[cohort["sample"] == m["sample"]]
        binders = sub[sub["binder_class"] != "Non-binder"]
        clonal = int((binders["clonality"] == "Clonal").sum())
        subclonal = int((binders["clonality"] == "Subclonal").sum())
        denom = clonal + subclonal
        rows.append({
            "sample": m["sample"], "patient": m["patient"],
            "timepoint": m["timepoint"],
            "total_peptides": int(len(sub)),
            "strong_binders": int((sub["binder_class"] == "Strong").sum()),
            "weak_binders": int((sub["binder_class"] == "Weak").sum()),
            "non_binders": int((sub["binder_class"] == "Non-binder").sum()),
            "clonal_binders": clonal, "subclonal_binders": subclonal,
            "fraction_clonal": round(clonal / denom, 4) if denom else np.nan,
        })
    return pd.DataFrame(rows)


def module5_neoepitope_nmd_summary(cohort: pd.DataFrame,
                                   variant_tbl: pd.DataFrame) -> pd.DataFrame:
    """Module 5 — per-sample neoepitopes split by their variant's NMD class.

    Inputs:  cohort table (augmented), and the per-variant NMD table (Module B)
             used to look up each neoepitope's variant-level NMD consensus.
    Output:  one row per sample with columns sample, patient, timepoint,
             total_neoepitopes_strong, nmd_sensitive_neo_strong,
             nmd_escape_neo_strong, total_neoepitopes_all, nmd_sensitive_neo_all,
             nmd_escape_neo_all, fraction_escape_strong, fraction_escape_all.
             A neoepitope is any peptide with binder_class != Non-binder; '_strong'
             restricts to Strong binders. fraction_escape =
             escape/(escape+sensitive); UNCERTAIN/UNKNOWN variants are excluded
             from both terms; NaN if the denominator is zero.
    """
    meta = _sample_meta(cohort)
    nmd_map = dict(zip(variant_tbl["vkey"], variant_tbl["nmd"]))
    df = cohort.copy()
    df["vkey"] = _variant_keys(df)
    df["vcons"] = df["vkey"].map(nmd_map)
    neo = df[df["binder_class"] != "Non-binder"]

    rows: List[Dict[str, object]] = []
    for _, m in meta.iterrows():
        sub = neo[neo["sample"] == m["sample"]]
        strong = sub[sub["binder_class"] == "Strong"]
        s_sens = int((strong["vcons"] == "SENSITIVE").sum())
        s_esc = int((strong["vcons"] == "INSENSITIVE").sum())
        a_sens = int((sub["vcons"] == "SENSITIVE").sum())
        a_esc = int((sub["vcons"] == "INSENSITIVE").sum())
        s_denom, a_denom = s_sens + s_esc, a_sens + a_esc
        rows.append({
            "sample": m["sample"], "patient": m["patient"],
            "timepoint": m["timepoint"],
            "total_neoepitopes_strong": int(len(strong)),
            "nmd_sensitive_neo_strong": s_sens,
            "nmd_escape_neo_strong": s_esc,
            "total_neoepitopes_all": int(len(sub)),
            "nmd_sensitive_neo_all": a_sens,
            "nmd_escape_neo_all": a_esc,
            "fraction_escape_strong": round(s_esc / s_denom, 4) if s_denom else np.nan,
            "fraction_escape_all": round(a_esc / a_denom, 4) if a_denom else np.nan,
        })
    return pd.DataFrame(rows)


# ═════════════════════════════════════════════════════════════════════════════
# PAIRED STATISTICS (Stage 3 — Wilcoxon signed-rank + BH-FDR)
# ═════════════════════════════════════════════════════════════════════════════

# (module, metric) pairs tested in paired_stats_stage3.tsv.
_STAGE3_PAIRED_SPEC: List[Tuple[str, str]] = [
    ("3", "nmd_fraction"), ("3", "nmd_sensitive"), ("3", "nmd_insensitive"),
    ("4", "total_peptides"), ("4", "strong_binders"),
    ("4", "clonal_binders"), ("4", "fraction_clonal"),
    ("5", "fraction_escape_strong"), ("5", "fraction_escape_all"),
    ("5", "nmd_sensitive_neo_strong"), ("5", "nmd_escape_neo_strong"),
]


def _bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR; NaN p-values are excluded and stay NaN."""
    out = np.full(pvals.shape, np.nan)
    mask = ~np.isnan(pvals)
    p = pvals[mask]
    n = p.size
    if n == 0:
        return out
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[mask] = adj
    return out


def compute_paired_stats_stage3(wide: pd.DataFrame) -> pd.DataFrame:
    """Paired Wilcoxon tests for the Module 3/4/5 metrics with BH-FDR.

    Input:   ``wide`` per-sample frame (sample, patient, timepoint + all module
             metric columns).
    Output:  DataFrame with columns module, metric, n_pairs, median_primary,
             median_recurrent, median_delta, mean_delta, wilcoxon_stat, p_value,
             p_value_adj_bh, note. Pairs require both timepoints present for a
             patient; pairs with a NaN metric value are dropped (n_pairs reflects
             usable pairs). BH-FDR is applied across all rows of this table.
    """
    prim = wide[wide.timepoint == "primary"].set_index("patient")
    recu = wide[wide.timepoint == "recurrent"].set_index("patient")
    paired = sorted(set(prim.index) & set(recu.index), key=int)

    rows: List[Dict[str, object]] = []
    for module, metric in _STAGE3_PAIRED_SPEC:
        if metric not in wide.columns:
            rows.append({"module": module, "metric": metric, "n_pairs": 0,
                         "median_primary": np.nan, "median_recurrent": np.nan,
                         "median_delta": np.nan, "mean_delta": np.nan,
                         "wilcoxon_stat": np.nan, "p_value": np.nan,
                         "p_value_adj_bh": np.nan, "note": "metric absent"})
            continue
        p = prim.loc[paired, metric].astype(float).to_numpy()
        r = recu.loc[paired, metric].astype(float).to_numpy()
        mask = ~(np.isnan(p) | np.isnan(r))
        pv, rv = p[mask], r[mask]
        delta = rv - pv
        n = len(delta)
        dropped = int(len(paired) - n)

        stat = p_value = np.nan
        note = f"{dropped} pair(s) dropped (undefined metric)" if dropped else ""
        if n < 6:
            note = (note + "; " if note else "") + "fewer than 6 pairs"
        elif np.allclose(delta, delta[0]):
            note = (note + "; " if note else "") + "zero variance in deltas"
        else:
            try:
                res = wilcoxon(rv, pv)
                stat, p_value = float(res.statistic), float(res.pvalue)
            except ValueError as exc:
                note = (note + "; " if note else "") + f"wilcoxon undefined ({exc})"

        rows.append({
            "module": module, "metric": metric, "n_pairs": n,
            "median_primary": round(float(np.median(pv)), 4) if n else np.nan,
            "median_recurrent": round(float(np.median(rv)), 4) if n else np.nan,
            "median_delta": round(float(np.median(delta)), 4) if n else np.nan,
            "mean_delta": round(float(np.mean(delta)), 4) if n else np.nan,
            "wilcoxon_stat": stat, "p_value": p_value,
            "p_value_adj_bh": np.nan, "note": note,
        })
    stats = pd.DataFrame(rows)
    stats["p_value_adj_bh"] = _bh_adjust(stats["p_value"].to_numpy(dtype=float))
    return stats


# ═════════════════════════════════════════════════════════════════════════════
# MODULE PLOTS (saved as 300 dpi PNGs to out_dir; also embedded in the report)
# ═════════════════════════════════════════════════════════════════════════════

# Palette shared with Stage 1.
COL_SENS = "#D85A30"   # NMD-sensitive
COL_ESC = "#378ADD"    # NMD-escape (insensitive)
COL_CLONAL = "#377AB8"
COL_SUBCLONAL = "#D85A30"


def _save_b64(fig: plt.Figure, out_dir: Path, name: str) -> str:
    """Save ``fig`` at 300 dpi to out_dir/name and return a base64 data URI."""
    path = out_dir / name
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()


def _paired_arrays(wide: pd.DataFrame, metric: str
                   ) -> Tuple[List[str], List[float], List[float]]:
    """Return (paired_patients, primary_values, recurrent_values) for ``metric``."""
    prim = wide[wide.timepoint == "primary"].set_index("patient")[metric]
    recu = wide[wide.timepoint == "recurrent"].set_index("patient")[metric]
    patients = sorted(set(prim.index) & set(recu.index), key=int)
    return (patients, [float(prim[p]) for p in patients],
            [float(recu[p]) for p in patients])


def _p_annotation(stats: pd.DataFrame, metric: str) -> str:
    """Format 'p = … (Δ̃ = …)' annotation for ``metric`` from the stats table."""
    row = stats[stats.metric == metric]
    if row.empty:
        return ""
    r = row.iloc[0]
    parts = []
    if pd.notna(r["p_value"]):
        p = float(r["p_value"])
        parts.append(f"p = {p:.1e}" if p < 1e-3 else f"p = {p:.3f}")
    else:
        parts.append("p = n/a")
    if pd.notna(r["median_delta"]):
        parts.append(f"median Δ = {float(r['median_delta']):.3g}")
    return "\n".join(parts)


def _connected_dots(ax: plt.Axes, prim: List[float], recu: List[float],
                    ylabel: str, title: str) -> None:
    """Draw a primary→recurrent connected dot plot onto ``ax``."""
    for yp, yr in zip(prim, recu):
        ax.plot([0, 1], [yp, yr], color="#cccccc", lw=0.9, zorder=1)
    ax.scatter([0] * len(prim), prim, color=COL_T, s=42, zorder=3)
    ax.scatter([1] * len(recu), recu, color=COL_M, s=42, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Primary", "Recurrent"])
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)


def plot_paired_metric(wide: pd.DataFrame, stats: pd.DataFrame, metric: str,
                       ylabel: str, title: str, out_dir: Path, name: str,
                       annotate: bool = True) -> str:
    """Connected dot plot for one paired metric (optionally p-value annotated)."""
    patients, prim, recu = _paired_arrays(wide, metric)
    fig, ax = plt.subplots(figsize=(6, 7))
    _connected_dots(ax, prim, recu, ylabel, title)
    if annotate:
        ann = _p_annotation(stats, metric)
        if ann:
            ax.text(0.97, 0.97, ann, transform=ax.transAxes, ha="right", va="top",
                    fontsize=11, bbox=dict(boxstyle="round", fc="white", ec="#cccccc"))
    return _save_b64(fig, out_dir, name)


def plot_waterfall(wide: pd.DataFrame, metric: str, ylabel: str, title: str,
                   out_dir: Path, name: str) -> str:
    """Waterfall of (recurrent − primary) per patient, sorted by delta."""
    patients, prim, recu = _paired_arrays(wide, metric)
    deltas = sorted([(p, r - t) for p, t, r in zip(patients, prim, recu)],
                    key=lambda x: x[1])
    labels = [d[0] for d in deltas]
    values = [d[1] for d in deltas]
    colors = [COL_SENS if v > 0 else COL_ESC for v in values]
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.bar(range(len(values)), values, color=colors, width=0.7)
    ax.axhline(0, color="#333333", lw=0.8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_xlabel("Patient (sorted by Δ)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    return _save_b64(fig, out_dir, name)


def plot_stacked_two(wide: pd.DataFrame, lower: str, upper: str,
                     lower_label: str, upper_label: str, lower_color: str,
                     upper_color: str, ylabel: str, title: str,
                     out_dir: Path, name: str) -> str:
    """Two-category stacked bar per sample (samples ordered patient↑, T before M)."""
    df = wide.sort_values(
        ["patient", "timepoint"],
        key=lambda s: s.map(int) if s.name == "patient"
        else s.map({"primary": 0, "recurrent": 1}))
    x = range(len(df))
    lo = df[lower].fillna(0).to_numpy()
    up = df[upper].fillna(0).to_numpy()
    fig, ax = plt.subplots(figsize=(16, 6))
    ax.bar(x, lo, color=lower_color, label=lower_label)
    ax.bar(x, up, bottom=lo, color=upper_color, label=upper_label)
    ax.set_xticks(list(x))
    ax.set_xticklabels(df["sample"], rotation=90, fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12)
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    return _save_b64(fig, out_dir, name)


def plot_scatter(wide: pd.DataFrame, xcol: str, ycol: str, xlabel: str,
                 ylabel: str, title: str, out_dir: Path, name: str,
                 logx: bool = False, spearman: bool = True) -> str:
    """Scatter of ycol vs xcol coloured by timepoint, optional Spearman ρ."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for tp, color in (("primary", COL_T), ("recurrent", COL_M)):
        sub = wide[wide.timepoint == tp]
        ax.scatter(sub[xcol], sub[ycol], color=color, s=42, alpha=0.8,
                   label=tp.capitalize())
    if spearman:
        x = pd.to_numeric(wide[xcol], errors="coerce").to_numpy()
        y = pd.to_numeric(wide[ycol], errors="coerce").to_numpy()
        ok = ~(np.isnan(x) | np.isnan(y))
        if logx:
            ok &= x > 0
        if ok.sum() > 2:
            rho, p = spearmanr(x[ok], y[ok])
            ann = (f"Spearman ρ = {rho:.2f}\n"
                   + (f"p = {p:.1e}" if p < 1e-3 else f"p = {p:.3f}"))
            ax.text(0.97, 0.97, ann, transform=ax.transAxes, ha="right", va="top",
                    fontsize=11, bbox=dict(boxstyle="round", fc="white", ec="#cccccc"))
    if logx:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12)
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    return _save_b64(fig, out_dir, name)


def plot_paired_two_panel(wide: pd.DataFrame, m1: str, m2: str, lab1: str,
                          lab2: str, title: str, out_dir: Path, name: str) -> str:
    """Two connected dot plots side by side (e.g. total vs strong binders)."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.5))
    for ax, metric, lab in ((axes[0], m1, lab1), (axes[1], m2, lab2)):
        _, prim, recu = _paired_arrays(wide, metric)
        _connected_dots(ax, prim, recu, lab, lab)
    fig.suptitle(title, fontsize=13)
    return _save_b64(fig, out_dir, name)


def make_module_plots(wide: pd.DataFrame, stats: pd.DataFrame, out_dir: Path,
                      have_stage1: bool) -> Dict[str, str]:
    """Generate the 12 Module 3/4/5 plots; return {key: base64 data URI}.

    The m4 TMB-vs-neoepitope scatter is skipped (empty URI) when Stage 1 data
    (tmb_per_mb) is unavailable.
    """
    print("[STEP] Generating Module 3/4/5 plots ...")
    b: Dict[str, str] = {}

    # Module 3
    b["m3_paired"] = plot_paired_metric(
        wide, stats, "nmd_fraction", "NMD fraction",
        "M3 — Paired NMD-sensitive fraction (FS)", out_dir,
        "m3_paired_nmd_fraction.png")
    b["m3_delta"] = plot_waterfall(
        wide, "nmd_fraction", "Δ NMD fraction (recurrent − primary)",
        "M3 — Change in NMD-sensitive fraction", out_dir,
        "m3_delta_nmd_fraction.png")
    b["m3_stacked"] = plot_stacked_two(
        wide, "nmd_sensitive", "nmd_insensitive", "NMD-sensitive",
        "NMD-escape (insensitive)", COL_SENS, COL_ESC, "FS variant count",
        "M3 — NMD-sensitive vs escape per sample", out_dir,
        "m3_stacked_sensitive_escape.png")
    b["m3_scatter"] = plot_scatter(
        wide, "truncating_in_pvacseq", "nmd_sensitive",
        "Truncating variants in pVACseq", "NMD-sensitive variants",
        "M3 — Truncating burden vs NMD-sensitive count", out_dir,
        "m3_scatter_trunc_vs_nmdsens.png")

    # Module 4
    b["m4_paired"] = plot_paired_two_panel(
        wide, "total_peptides", "strong_binders", "Total binders",
        "Strong binders", "M4 — Paired neoepitope burden", out_dir,
        "m4_paired_neoepitope_burden.png")
    b["m4_stacked"] = plot_stacked_two(
        wide, "clonal_binders", "subclonal_binders", "Clonal binders",
        "Subclonal binders", COL_CLONAL, COL_SUBCLONAL, "Binder count",
        "M4 — Binder clonality per sample", out_dir, "m4_stacked_clonality.png")
    if have_stage1 and "tmb_per_mb" in wide.columns:
        b["m4_scatter"] = plot_scatter(
            wide, "tmb_per_mb", "strong_binders", "TMB per Mb (log scale)",
            "Strong binders", "M4 — TMB vs strong-binder burden", out_dir,
            "m4_scatter_tmb_vs_neo.png", logx=True)
    else:
        print("  [WARN] Stage 1 TMB unavailable — skipping m4_scatter_tmb_vs_neo.png")
        b["m4_scatter"] = ""
    b["m4_delta"] = plot_waterfall(
        wide, "strong_binders", "Δ strong binders (recurrent − primary)",
        "M4 — Change in strong-binder burden", out_dir,
        "m4_delta_neoepitopes.png")

    # Module 5
    b["m5_paired"] = plot_paired_metric(
        wide, stats, "fraction_escape_strong", "Fraction NMD-escape (strong)",
        "M5 — Paired NMD-escape fraction (strong binders)", out_dir,
        "m5_paired_fraction_escape.png")
    b["m5_stacked"] = plot_stacked_two(
        wide, "nmd_sensitive_neo_strong", "nmd_escape_neo_strong",
        "NMD-sensitive neo", "NMD-escape neo", COL_SENS, COL_ESC,
        "Strong neoepitope count",
        "M5 — Sensitive vs escape neoepitopes per sample (strong)", out_dir,
        "m5_stacked_sensitive_escape_neo.png")
    b["m5_delta"] = plot_waterfall(
        wide, "fraction_escape_strong",
        "Δ NMD-escape fraction (recurrent − primary)",
        "M5 — Change in NMD-escape fraction (strong)", out_dir,
        "m5_delta_fraction_escape.png")
    b["m5_paired_count"] = plot_paired_metric(
        wide, stats, "nmd_sensitive_neo_strong", "NMD-sensitive neoepitopes",
        "M5 — Paired NMD-sensitive neoepitope count (strong)", out_dir,
        "m5_paired_nmd_sensitive_neo.png")
    return b


# ═════════════════════════════════════════════════════════════════════════════
# PER-PAIR VARIANT OVERLAP WIDGET DATA
# ═════════════════════════════════════════════════════════════════════════════

def build_overlap_data(overlap_df: Optional[pd.DataFrame],
                       variant_tbl: pd.DataFrame,
                       stage1: Optional[pd.DataFrame]) -> Dict[str, dict]:
    """Build the per-pair overlap widget JSON from Stage 1's overlap table.

    Inputs:  ``overlap_df`` Stage 1 paired_variant_overlap.tsv (or None); the
             per-variant NMD table; and the optional Stage 1 summary. Patients
             are kept only if both timepoints carry ≥1 truncating variant. The
             truncating set is taken from Stage 1's truncating_total (whole-exome
             FS+SG+splice+start_lost) when available, else falls back to FS
             variants that reached pVACseq.
    Output:  {patient: {t_total, m_total, shared, t_private, m_private,
             shared_pct_of_t, shared_pct_of_m}}. Empty dict if overlap data is
             unavailable.
    """
    if overlap_df is None or overlap_df.empty:
        return {}
    if stage1 is not None and "truncating_total" in stage1.columns:
        s = stage1.copy()
        s["patient"] = s["sample"].astype(str).str.split("_").str[0]
        s["tp"] = s["sample"].astype(str).str.split("_").str[1]
        trunc = pd.to_numeric(s["truncating_total"], errors="coerce").fillna(0)
        has_t = set(s[(s.tp == "T") & (trunc > 0)]["patient"])
        has_m = set(s[(s.tp == "M") & (trunc > 0)]["patient"])
    else:
        fs = variant_tbl[variant_tbl["Variant Type"] == "FS"]
        has_t = set(fs[fs.timepoint == "primary"]["patient"].astype(str))
        has_m = set(fs[fs.timepoint == "recurrent"]["patient"].astype(str))
    qualifying = has_t & has_m
    out: Dict[str, dict] = {}
    for _, r in overlap_df.iterrows():
        pid = str(int(r["patient"])) if not isinstance(r["patient"], str) else str(r["patient"])
        if qualifying and pid not in qualifying:
            continue
        out[pid] = {
            "t_total": int(r["t_total"]), "m_total": int(r["m_total"]),
            "shared": int(r["shared"]), "t_private": int(r["t_private"]),
            "m_private": int(r["m_private"]),
            "shared_pct_of_t": float(r["shared_pct_of_t"]),
            "shared_pct_of_m": float(r["shared_pct_of_m"]),
        }
    return out


# ═════════════════════════════════════════════════════════════════════════════
# PER-SAMPLE WIDGET DATA (JSON-serialisable for embedding in HTML)
# ═════════════════════════════════════════════════════════════════════════════

# Per-algorithm IC50 columns (pVACseq output column names)
_ALGO_IC50_COLS = {
    "MHCflurry":    "MHCflurry MT IC50 Score",
    "MHCnuggetsI":  "MHCnuggetsI MT IC50 Score",
    "NetMHC":       "NetMHC MT IC50 Score",
    "NetMHCpan":    "NetMHCpan MT IC50 Score",
    "PickPocket":   "PickPocket MT IC50 Score",
    "SMM":          "SMM MT IC50 Score",
    "SMMPMBEC":     "SMMPMBEC MT IC50 Score",
}


def _num(v):
    """Coerce a TSV cell to float-or-None; tolerates NaN, empty, and 'NA'."""
    try:
        if pd.isna(v): return None
        x = float(v)
        return None if pd.isna(x) else round(x, 3)
    except (TypeError, ValueError):
        return None


def _candidate_dict(r) -> dict:
    """Convert one row of the per-sample TSV into a JS-friendly dict for the
    interactive widget. Includes per-algorithm IC50 scores so the candidate
    dropdown can compare what each binding predictor said."""
    methods = {label: _num(r.get(col)) for label, col in _ALGO_IC50_COLS.items()}
    return {
        "gene":   str(r.get("Gene Name", "")),
        "pep":    str(r.get("MT Epitope Seq", "")),
        "hla":    str(r.get("hla_allele", "")),
        "vtype":  str(r.get("Variant Type", "")),
        "ic50":   _num(r.get("best_ic50")),
        "median": _num(r.get("median_ic50")) if "median_ic50" in r else _num(r.get("Median MT IC50 Score")),
        "method": str(r.get("best_ic50_method", "")),
        "methods": methods,
        "nmd":    str(r.get("nmd_consensus", "")),
        "conf":   str(r.get("nmd_confidence", "")),
        "tier":   str(r.get("priority_tier", "")),
        "rule":   str(r.get("nmd_rule_explanation", "")),
    }



def build_sample_data(fs_cohort: pd.DataFrame) -> dict:
    """Per-sample data for the interactive drill-down widget.
    Filtered to FS only (NMD-actionable). Empty samples retained as empty dicts."""
    out = {}
    if fs_cohort.empty:
        return out
    for sample, sub in fs_cohort.groupby("sample"):
        sub = sub.sort_values("best_ic50", na_position="last")
        out[sample] = {
            "patient":   sub["patient"].iloc[0],
            "timepoint": sub["timepoint"].iloc[0],
            "n_total":   int(len(sub)),
            "tier1":     int((sub["priority_tier"] == "TIER1").sum()),
            "tier2":     int((sub["priority_tier"] == "TIER2").sum()),
            "tier3":     int((sub["priority_tier"] == "TIER3_control").sum()),
            "unclass":   int((sub["priority_tier"] == "UNCLASSIFIED").sum()),
            "sensitive": int((sub["nmd_consensus"] == "SENSITIVE").sum()),
            "insensitive": int((sub["nmd_consensus"] == "INSENSITIVE").sum()),
            "candidates": [_candidate_dict(r) for _, r in sub.iterrows()]
        }
    return out


def build_patient_data(fs_cohort: pd.DataFrame) -> dict:
    """Per-patient data for the interactive drill-down widget.
    Each patient has T (primary) and M (recurrent) panels - either may be None
    if that timepoint has no FS candidates. FS only (NMD-actionable)."""
    out = {}
    if fs_cohort.empty:
        return out

    def _panel(sub):
        sub = sub.sort_values("best_ic50", na_position="last")
        return {
            "n_total":     int(len(sub)),
            "tier1":       int((sub["priority_tier"] == "TIER1").sum()),
            "tier2":       int((sub["priority_tier"] == "TIER2").sum()),
            "tier3":       int((sub["priority_tier"] == "TIER3_control").sum()),
            "unclass":     int((sub["priority_tier"] == "UNCLASSIFIED").sum()),
            "sensitive":   int((sub["nmd_consensus"] == "SENSITIVE").sum()),
            "insensitive": int((sub["nmd_consensus"] == "INSENSITIVE").sum()),
            "candidates":  [_candidate_dict(r) for _, r in sub.iterrows()],
        }

    for patient, sub in fs_cohort.groupby("patient"):
        t_sub = sub[sub["timepoint"] == "primary"]
        m_sub = sub[sub["timepoint"] == "recurrent"]
        out[str(patient)] = {
            "patient": str(patient),
            "T": _panel(t_sub) if len(t_sub) else None,
            "M": _panel(m_sub) if len(m_sub) else None,
        }
    return out


# ═════════════════════════════════════════════════════════════════════════════
# HTML REPORT
# ═════════════════════════════════════════════════════════════════════════════

CSS = """
*{box-sizing:border-box} body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;font-size:14px;line-height:1.6;color:#1a1a1a;background:#fff;margin:0;padding:0}
.page{max-width:1200px;margin:0 auto;padding:40px 32px}
h1{font-size:24px;font-weight:500;margin-bottom:4px}
h2{font-size:17px;font-weight:500;border-bottom:1px solid #e0e0e0;padding-bottom:6px;margin:36px 0 8px}
h3{font-size:14px;font-weight:500;margin:16px 0 6px;color:#444}
p.sub{color:#666;font-size:13px;margin:0 0 16px;line-height:1.5}
table{border-collapse:collapse;width:100%;font-size:12px;margin-top:8px}
th{background:#f4f4f4;padding:6px 10px;text-align:left;font-weight:500;border-bottom:2px solid #ddd}
td{padding:6px 10px;border-bottom:1px solid #eee}
tr:hover td{background:#fafafa}
.card{background:#f8f8f8;border-radius:8px;padding:14px 18px;display:inline-block;min-width:140px;margin:4px;vertical-align:top}
.cv{font-size:22px;font-weight:500} .cl{font-size:11px;color:#777}
.cs{font-size:11px;color:#999;margin-top:2px}
.note{background:#fffbe6;border-left:3px solid #f0c040;padding:10px 14px;margin:12px 0;font-size:13px;color:#555}
.intro-box{background:#f4f8fc;border-left:3px solid #378ADD;padding:14px 18px;margin:12px 0;font-size:13px;color:#444;line-height:1.6}
.tier-block{background:#fafafa;padding:10px 14px;margin:8px 0;border-radius:4px;font-size:13px}
.tier-1-h{color:#B8860B;font-weight:500}
.tier-2-h{color:#C57A33;font-weight:500}
.tier-3-h{color:#377AB8;font-weight:500}
.widget{background:#fafafa;border:1px solid #e0e0e0;border-radius:8px;padding:18px 22px;margin:16px 0}
.widget select{font-size:14px;padding:6px 10px;border:1px solid #ccc;border-radius:4px;background:#fff;font-family:inherit}
.widget label{font-weight:500;margin-right:10px}
.widget .meta{font-size:12px;color:#777;margin-top:4px}
"""


def _card(label, val, sub=""):
    s = f'<div class="cs">{sub}</div>' if sub else ""
    return (f'<div class="card"><div class="cl">{label}</div>'
            f'<div class="cv">{val}</div>{s}</div>')


def _tbl(df, cols, empty="No candidates."):
    if df is None or df.empty:
        return f"<p style='color:#888;font-style:italic;'>{empty}</p>"
    th = "".join(f"<th>{c}</th>" for c in cols)
    rows = []
    for _, r in df.iterrows():
        rows.append("<tr>" + "".join(f"<td>{r.get(c,'')}</td>" for c in cols) + "</tr>")
    return f"<table><thead><tr>{th}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def _img(b64, caption):
    if not b64:
        return ""
    return (f'<div style="margin:16px 0;"><img src="{b64}" '
            f'style="max-width:90%;display:block;margin:0 auto;border-radius:4px;'
            f'box-shadow:0 1px 4px rgba(0,0,0,.12);">'
            f'<p style="font-size:11px;color:#888;margin-top:6px;text-align:center;">{caption}</p></div>')


METHODS_BLOCK = """
<h2>NMD scoring methods</h2>
<h3>Method 1 — VEP NMD plugin</h3>
<p class="sub">
  The VEP NMD plugin (Ensembl v105+) annotates truncating variants with
  <code>NMD_escaping_variant</code> when the PTC is predicted to escape NMD. An empty
  NMD field on a truncating variant means NMD is triggered (SENSITIVE).
</p>

<h3>Method 2 — Lindeboom rules (Nat Genet 2019)</h3>
<p class="sub">Applied to frameshift and stop-gained variants. Rules in priority order:</p>
<ul style="font-size:13px;color:#444;margin:0 0 12px 20px;line-height:1.8">
  <li><strong>Rule 4 — Start-proximal PTC (&lt;150 nt):</strong> Pioneer round completes before NMD surveillance → INSENSITIVE</li>
  <li><strong>Rule 1 — Last exon:</strong> No downstream EJC to trigger NMD → INSENSITIVE</li>
  <li><strong>Rule 3 — Long exon (&gt;407 nt):</strong> EJC too far downstream → INSENSITIVE</li>
  <li><strong>Rule 2 — 55 nt boundary:</strong> PTC &gt;55 nt upstream of last EJC → SENSITIVE (canonical NMD)</li>
</ul>

<h3>Ensemble confidence (0–3)</h3>
<p class="sub">3 = both methods agree; 2 = single method available; 1 = methods disagree; 0 = no data.</p>

<h2>Priority tiers</h2>
<div class="tier-block"><span class="tier-1-h">TIER 1</span> — NMD-sensitive + IC50 &lt; 50 nM. Primary therapeutic targets — silenced by NMD, exposed by NMD inhibition.</div>
<div class="tier-block"><span class="tier-2-h">TIER 2</span> — NMD-sensitive + IC50 50–500 nM. Moderate binders, potentially relevant after NMD inhibition.</div>
<div class="tier-block"><span class="tier-3-h">TIER 3</span> — NMD-insensitive + IC50 &lt; 500 nM. Already expressed — controls for immune response without NMD inhibition.</div>
"""


WIDGET_HTML = """
<h2>Per-patient drill-down (interactive)</h2>
<p class="sub">
  Select a patient to see their NMD-actionable (frameshift) candidates in both
  primary (T) and recurrent (M) timepoints side by side. Then pick a candidate
  in either panel to see how the 7 binding-prediction algorithms compare.
  NetMHCpanEL is excluded because it reports an eluted-ligand probability rather than IC50.
</p>
<div class="widget">
  <div style="display:flex;gap:14px;align-items:center;margin-bottom:14px;">
    <label for="patient-select"><strong>Patient:</strong></label>
    <select id="patient-select" style="min-width:240px;"></select>
  </div>

  <div style="display:grid;grid-template-columns:1fr 1fr;gap:18px;">
    <div class="tp-panel" style="border:1px solid #e0e0e0;border-radius:6px;padding:14px;min-width:0;">
      <h3 style="margin-top:0;margin-bottom:8px;color:#377AB8;border-bottom:2px solid #377AB8;padding-bottom:4px;">Primary (T)</h3>
      <div class="meta" id="meta-T" style="font-size:12px;"></div>
      <div id="bars-T" style="margin:8px 0;"></div>
      <div style="margin-top:8px;">
        <label for="cand-T" style="font-size:12px;">Candidate:</label><br>
        <select id="cand-T" style="width:100%;margin-top:4px;"></select>
      </div>
      <div class="meta" id="cmeta-T" style="margin-top:10px;font-size:12px;"></div>
      <div id="methods-T" style="margin-top:8px;"></div>
    </div>
    <div class="tp-panel" style="border:1px solid #e0e0e0;border-radius:6px;padding:14px;min-width:0;">
      <h3 style="margin-top:0;margin-bottom:8px;color:#D85A30;border-bottom:2px solid #D85A30;padding-bottom:4px;">Recurrent (M)</h3>
      <div class="meta" id="meta-M" style="font-size:12px;"></div>
      <div id="bars-M" style="margin:8px 0;"></div>
      <div style="margin-top:8px;">
        <label for="cand-M" style="font-size:12px;">Candidate:</label><br>
        <select id="cand-M" style="width:100%;margin-top:4px;"></select>
      </div>
      <div class="meta" id="cmeta-M" style="margin-top:10px;font-size:12px;"></div>
      <div id="methods-M" style="margin-top:8px;"></div>
    </div>
  </div>

  <h3 style="margin-top:24px;">All FS candidates for this patient (both timepoints)</h3>
  <div id="patient-table" style="margin-top:6px;"></div>
</div>
<script>
const PATIENT_DATA = __PATIENT_DATA__;
const PATIENT_KEYS = Object.keys(PATIENT_DATA).sort((a,b) => parseInt(a) - parseInt(b));

function svgBarV(items, w=380, h=180, yLabel='') {
  const max = Math.max(1, ...items.map(x => x.value));
  const padL = 50, padR = 12, padT = 16, padB = 36;
  const innerW = w - padL - padR, innerH = h - padT - padB;
  const barW = innerW / items.length - 10;
  let svg = `<svg viewBox="0 0 ${w} ${h}" xmlns="http://www.w3.org/2000/svg" style="font-family:inherit;font-size:11px;max-width:100%;height:auto;">`;
  svg += `<line x1="${padL}" y1="${padT}" x2="${padL}" y2="${h-padB}" stroke="#bbb"/>`;
  svg += `<line x1="${padL}" y1="${h-padB}" x2="${w-padR}" y2="${h-padB}" stroke="#bbb"/>`;
  for (let i=0; i<=4; i++) {
    const y = h - padB - (i/4)*innerH;
    const v = Math.round(max * i / 4 * 10) / 10;
    svg += `<line x1="${padL-3}" y1="${y}" x2="${padL}" y2="${y}" stroke="#bbb"/>`;
    svg += `<text x="${padL-6}" y="${y+3}" text-anchor="end" fill="#666">${v}</text>`;
  }
  items.forEach((it, i) => {
    const bh = (it.value / max) * innerH;
    const x = padL + i * (barW + 10) + 5;
    const y = h - padB - bh;
    svg += `<rect x="${x}" y="${y}" width="${barW}" height="${bh}" fill="${it.color}" rx="2"/>`;
    if (it.value > 0) {
      svg += `<text x="${x + barW/2}" y="${y - 4}" text-anchor="middle" fill="#333">${it.value}</text>`;
    }
    svg += `<text x="${x + barW/2}" y="${h - padB + 14}" text-anchor="middle" fill="#444">${it.label}</text>`;
  });
  svg += '</svg>';
  return svg;
}

function svgBarH_log(items, w=480, h=240, xLabel='IC50 (nM, log scale)') {
  const padL = 110, padR = 50, padT = 16, padB = 42;
  const innerW = w - padL - padR;
  const innerH = h - padT - padB;
  const barH = Math.max(8, innerH / items.length - 6);

  const vals = items.map(x => x.value).filter(v => v != null && v > 0);
  if (vals.length === 0) {
    return `<p style="color:#888;font-style:italic;">No IC50 values to display.</p>`;
  }
  const dataMin = Math.max(0.01, Math.min(...vals));
  const dataMax = Math.max(...vals);
  const lmin = Math.log10(dataMin * 0.7);
  const lmax = Math.log10(Math.max(dataMax * 1.3, 1000));
  const xScale = v => padL + (Math.log10(Math.max(0.01, v)) - lmin) / (lmax - lmin) * innerW;

  let svg = `<svg viewBox="0 0 ${w} ${h}" xmlns="http://www.w3.org/2000/svg" style="font-family:inherit;font-size:10px;max-width:100%;height:auto;">`;
  svg += `<line x1="${padL}" y1="${padT}" x2="${padL}" y2="${h-padB}" stroke="#bbb"/>`;
  svg += `<line x1="${padL}" y1="${h-padB}" x2="${w-padR}" y2="${h-padB}" stroke="#bbb"/>`;

  for (let log = Math.ceil(lmin); log <= Math.floor(lmax); log++) {
    const v = Math.pow(10, log);
    const x = xScale(v);
    svg += `<line x1="${x}" y1="${h-padB}" x2="${x}" y2="${h-padB+3}" stroke="#bbb"/>`;
    svg += `<text x="${x}" y="${h-padB+14}" text-anchor="middle" fill="#666">${v >= 1 ? v : v.toFixed(1)}</text>`;
  }
  svg += `<text x="${padL + innerW/2}" y="${h-4}" text-anchor="middle" fill="#666">${xLabel}</text>`;

  if (Math.pow(10, lmin) <= 50 && 50 <= Math.pow(10, lmax)) {
    const x = xScale(50);
    svg += `<line x1="${x}" y1="${padT}" x2="${x}" y2="${h-padB}" stroke="#555" stroke-dasharray="3 3" stroke-width="0.8"/>`;
    svg += `<text x="${x+3}" y="${padT+10}" font-size="9" fill="#555">50 nM</text>`;
  }
  if (Math.pow(10, lmin) <= 500 && 500 <= Math.pow(10, lmax)) {
    const x = xScale(500);
    svg += `<line x1="${x}" y1="${padT}" x2="${x}" y2="${h-padB}" stroke="#999" stroke-dasharray="2 4" stroke-width="0.8"/>`;
    svg += `<text x="${x+3}" y="${padT+10}" font-size="9" fill="#999">500 nM</text>`;
  }

  items.forEach((it, i) => {
    const yMid = padT + i * (innerH / items.length) + (innerH / items.length) / 2;
    const yTop = yMid - barH/2;
    if (it.value == null || it.value <= 0) {
      svg += `<text x="${padL - 6}" y="${yMid+3}" text-anchor="end" fill="#444">${it.label}</text>`;
      svg += `<text x="${padL + 4}" y="${yMid+3}" fill="#aaa" font-style="italic">no value</text>`;
      return;
    }
    const xEnd = xScale(it.value);
    svg += `<text x="${padL - 6}" y="${yMid+3}" text-anchor="end" fill="#444">${it.label}</text>`;
    svg += `<rect x="${padL}" y="${yTop}" width="${xEnd - padL}" height="${barH}" fill="${it.color}" rx="2"/>`;
    svg += `<text x="${xEnd + 4}" y="${yMid+3}" fill="#333">${it.value < 1 ? it.value.toFixed(2) : it.value.toFixed(1)}</text>`;
  });

  svg += '</svg>';
  return svg;
}

function escapeHTML(s) {
  return String(s).replace(/[&<>"']/g, m => ({
    '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'
  }[m]));
}

function tierColor(t) {
  return {TIER1:'#F0C040',TIER2:'#F5944D',TIER3_control:'#378ADD',UNCLASSIFIED:'#bbbbbb'}[t] || '#999';
}

function nmdColor(n) {
  return {SENSITIVE:'#D85A30',INSENSITIVE:'#378ADD',UNCERTAIN:'#aaa',UNKNOWN:'#ccc'}[n] || '#999';
}

function renderTable(rows) {
  if (!rows || !rows.length) return '<p style="color:#888;font-style:italic;">No frameshift candidates for this patient.</p>';
  let html = '<table><thead><tr>'
    + ['Timepoint','Gene','Peptide','HLA','best IC50','median','best method','NMD','tier','rule'].map(h => `<th>${h}</th>`).join('')
    + '</tr></thead><tbody>';
  for (const r of rows) {
    const tierBadge = `<span style="display:inline-block;background:${tierColor(r.tier)};color:#fff;padding:2px 6px;border-radius:3px;font-size:11px;">${r.tier}</span>`;
    const tpColor = r.tp === 'T' ? '#377AB8' : '#D85A30';
    html += '<tr>'
      + `<td><strong style="color:${tpColor};">${r.tp}</strong></td>`
      + `<td>${escapeHTML(r.gene)}</td>`
      + `<td><code>${escapeHTML(r.pep)}</code></td>`
      + `<td>${escapeHTML(r.hla)}</td>`
      + `<td>${r.ic50 == null ? '—' : r.ic50.toFixed(2)}</td>`
      + `<td>${r.median == null ? '—' : r.median.toFixed(1)}</td>`
      + `<td>${escapeHTML(r.method)}</td>`
      + `<td style="color:${nmdColor(r.nmd)};">${escapeHTML(r.nmd)}</td>`
      + `<td>${tierBadge}</td>`
      + `<td style="font-size:11px;color:#666;">${escapeHTML(r.rule)}</td>`
      + '</tr>';
  }
  html += '</tbody></table>';
  return html;
}

const METHOD_ORDER = ['MHCflurry','MHCnuggetsI','NetMHC','NetMHCpan','PickPocket','SMM','SMMPMBEC'];

function renderPanel(pid, tp) {
  const d = PATIENT_DATA[pid];
  const panel = d ? d[tp] : null;
  const tpName = tp === 'T' ? 'primary' : 'recurrent';
  const metaEl = document.getElementById('meta-' + tp);
  const barsEl = document.getElementById('bars-' + tp);
  const candSel = document.getElementById('cand-' + tp);
  const cmetaEl = document.getElementById('cmeta-' + tp);
  const methodsEl = document.getElementById('methods-' + tp);

  candSel.innerHTML = '';
  if (!panel) {
    metaEl.innerHTML = `<em style="color:#aaa;">No FS candidates in ${tpName} timepoint.</em>`;
    barsEl.innerHTML = '';
    cmetaEl.innerHTML = '';
    methodsEl.innerHTML = '';
    candSel.disabled = true;
    return;
  }
  candSel.disabled = false;
  metaEl.innerHTML = `<strong>${panel.n_total}</strong> FS · `
    + `<span style="color:#B8860B;">${panel.tier1} T1</span> · `
    + `<span style="color:#C57A33;">${panel.tier2} T2</span> · `
    + `<span style="color:#377AB8;">${panel.tier3} T3 ctrl</span> · `
    + `<span style="color:#999;">${panel.unclass} unclass</span>`;
  barsEl.innerHTML = svgBarV([
    {label:'TIER1', value:panel.tier1, color:'#F0C040'},
    {label:'TIER2', value:panel.tier2, color:'#F5944D'},
    {label:'TIER3', value:panel.tier3, color:'#378ADD'},
    {label:'Unclass', value:panel.unclass, color:'#bbbbbb'},
  ]);

  panel.candidates.forEach((c, i) => {
    const opt = document.createElement('option');
    opt.value = i;
    opt.textContent = `${c.gene} · ${c.hla} · ${c.tier} · ${c.ic50 != null ? c.ic50.toFixed(2)+' nM' : '—'}`;
    candSel.appendChild(opt);
  });
  if (panel.candidates.length) renderCandidate(pid, tp, 0);
}

function renderCandidate(pid, tp, idx) {
  const d = PATIENT_DATA[pid];
  const panel = d ? d[tp] : null;
  if (!panel) return;
  const c = panel.candidates[idx];
  if (!c) return;
  const cmetaEl = document.getElementById('cmeta-' + tp);
  const methodsEl = document.getElementById('methods-' + tp);
  const tierBadge = `<span style="display:inline-block;background:${tierColor(c.tier)};color:#fff;padding:2px 6px;border-radius:3px;">${c.tier}</span>`;
  cmetaEl.innerHTML = `<strong>${c.gene}</strong> · <code>${c.pep}</code> · ${c.hla} · `
    + `${c.vtype} · ${tierBadge}<br>`
    + `NMD: <span style="color:${nmdColor(c.nmd)};">${c.nmd}</span> (${c.conf}) · `
    + `best: ${c.method}@${c.ic50 != null ? c.ic50.toFixed(2) : '—'} nM · median ${c.median != null ? c.median.toFixed(1) : '—'} nM`;
  if (c.rule) cmetaEl.innerHTML += `<br><span style="font-size:11px;color:#666;">${escapeHTML(c.rule)}</span>`;

  const items = METHOD_ORDER.map(m => ({
    label: m,
    value: c.methods[m] != null ? c.methods[m] : null,
    color: m === c.method ? '#D85A30' : '#378ADD'
  }));
  methodsEl.innerHTML = svgBarH_log(items);
}

function renderPatient(pid) {
  const d = PATIENT_DATA[pid];
  if (!d) return;
  renderPanel(pid, 'T');
  renderPanel(pid, 'M');
  const allCands = [];
  if (d.T) for (const c of d.T.candidates) allCands.push(Object.assign({}, c, {tp:'T'}));
  if (d.M) for (const c of d.M.candidates) allCands.push(Object.assign({}, c, {tp:'M'}));
  document.getElementById('patient-table').innerHTML = renderTable(allCands);
}

(function init() {
  const sel = document.getElementById('patient-select');
  for (const k of PATIENT_KEYS) {
    const d = PATIENT_DATA[k];
    const nT = d.T ? d.T.n_total : 0;
    const nM = d.M ? d.M.n_total : 0;
    const opt = document.createElement('option');
    opt.value = k;
    opt.textContent = `Patient ${k}  (T=${nT} · M=${nM} FS)`;
    sel.appendChild(opt);
  }
  sel.addEventListener('change', e => renderPatient(e.target.value));
  document.getElementById('cand-T').addEventListener('change', e =>
    renderCandidate(sel.value, 'T', parseInt(e.target.value)));
  document.getElementById('cand-M').addEventListener('change', e =>
    renderCandidate(sel.value, 'M', parseInt(e.target.value)));
  if (PATIENT_KEYS.length) renderPatient(PATIENT_KEYS[0]);
})();
</script>
"""


OVERLAP_WIDGET_HTML = """
<h2>Per-pair variant overlap (interactive)</h2>
<p class="sub">
  Select a patient to compare their primary (T) and recurrent (M) PASS variant
  sets (Stage 1 <code>paired_variant_overlap.tsv</code>). Only patients with at
  least one truncating (frameshift) variant in <em>both</em> timepoints are shown.
</p>
<div class="widget">
  <div style="display:flex;gap:14px;align-items:center;margin-bottom:14px;">
    <label for="overlap-select"><strong>Patient:</strong></label>
    <select id="overlap-select" style="min-width:200px;"></select>
  </div>
  <div style="display:flex;gap:28px;flex-wrap:wrap;align-items:center;">
    <div id="overlap-venn"></div>
    <div id="overlap-box" class="meta" style="font-size:13px;"></div>
  </div>
</div>
<script>
const OVERLAP_DATA = __OVERLAP_DATA__;
const OVERLAP_KEYS = Object.keys(OVERLAP_DATA).sort((a,b) => parseInt(a) - parseInt(b));

function overlapVenn(d, w=420, h=240) {
  const cx1 = 150, cx2 = 270, cy = 120, r = 95;
  const COL_T = '#377AB8', COL_M = '#D85A30';
  let svg = `<svg viewBox="0 0 ${w} ${h}" xmlns="http://www.w3.org/2000/svg" style="font-family:inherit;font-size:12px;max-width:100%;height:auto;">`;
  svg += `<circle cx="${cx1}" cy="${cy}" r="${r}" fill="${COL_T}" fill-opacity="0.32" stroke="${COL_T}" stroke-width="1.5"/>`;
  svg += `<circle cx="${cx2}" cy="${cy}" r="${r}" fill="${COL_M}" fill-opacity="0.32" stroke="${COL_M}" stroke-width="1.5"/>`;
  const pctT = d.t_total ? (100*d.t_private/d.t_total).toFixed(0) : 0;
  const pctM = d.m_total ? (100*d.m_private/d.m_total).toFixed(0) : 0;
  const pctST = d.shared_pct_of_t.toFixed(0), pctSM = d.shared_pct_of_m.toFixed(0);
  svg += `<text x="${cx1-42}" y="${cy-4}" text-anchor="middle" fill="#1a4a7a" font-weight="600">${d.t_private}</text>`;
  svg += `<text x="${cx1-42}" y="${cy+12}" text-anchor="middle" fill="#1a4a7a" font-size="10">T-only (${pctT}%)</text>`;
  svg += `<text x="${cx2+42}" y="${cy-4}" text-anchor="middle" fill="#9a3a14" font-weight="600">${d.m_private}</text>`;
  svg += `<text x="${cx2+42}" y="${cy+12}" text-anchor="middle" fill="#9a3a14" font-size="10">M-only (${pctM}%)</text>`;
  svg += `<text x="${(cx1+cx2)/2}" y="${cy-4}" text-anchor="middle" fill="#222" font-weight="600">${d.shared}</text>`;
  svg += `<text x="${(cx1+cx2)/2}" y="${cy+12}" text-anchor="middle" fill="#222" font-size="10">shared</text>`;
  svg += `<text x="${cx1}" y="${cy-r-8}" text-anchor="middle" fill="${COL_T}" font-weight="600">Primary (T): ${d.t_total}</text>`;
  svg += `<text x="${cx2}" y="${cy+r+20}" text-anchor="middle" fill="${COL_M}" font-weight="600">Recurrent (M): ${d.m_total}</text>`;
  svg += `<text x="${(cx1+cx2)/2}" y="${h-6}" text-anchor="middle" fill="#888" font-size="10">shared = ${pctST}% of T · ${pctSM}% of M</text>`;
  svg += '</svg>';
  return svg;
}

function renderOverlap(pid) {
  const d = OVERLAP_DATA[pid];
  if (!d) return;
  document.getElementById('overlap-venn').innerHTML = overlapVenn(d);
  document.getElementById('overlap-box').innerHTML =
      '<table style="font-size:12px;"><tbody>'
    + `<tr><td>t_total</td><td style="text-align:right;"><strong>${d.t_total}</strong></td></tr>`
    + `<tr><td>m_total</td><td style="text-align:right;"><strong>${d.m_total}</strong></td></tr>`
    + `<tr><td>shared</td><td style="text-align:right;"><strong>${d.shared}</strong></td></tr>`
    + `<tr><td>t_private</td><td style="text-align:right;"><strong>${d.t_private}</strong></td></tr>`
    + `<tr><td>m_private</td><td style="text-align:right;"><strong>${d.m_private}</strong></td></tr>`
    + `<tr><td>shared_pct_of_t</td><td style="text-align:right;"><strong>${d.shared_pct_of_t}%</strong></td></tr>`
    + `<tr><td>shared_pct_of_m</td><td style="text-align:right;"><strong>${d.shared_pct_of_m}%</strong></td></tr>`
    + '</tbody></table>';
}

(function initOverlap() {
  const sel = document.getElementById('overlap-select');
  if (!OVERLAP_KEYS.length) {
    document.getElementById('overlap-box').innerHTML =
      '<em style="color:#aaa;">Overlap data unavailable (Stage 1 paired_variant_overlap.tsv not provided).</em>';
    return;
  }
  for (const k of OVERLAP_KEYS) {
    const opt = document.createElement('option');
    opt.value = k; opt.textContent = `Patient ${k}`;
    sel.appendChild(opt);
  }
  sel.addEventListener('change', e => renderOverlap(e.target.value));
  renderOverlap(OVERLAP_KEYS[0]);
})();
</script>
"""


def _stats_table_html(stats: pd.DataFrame, module: str) -> str:
    """Render the paired_stats rows for one module as a compact HTML table."""
    sub = stats[stats["module"] == module].copy()
    if sub.empty:
        return "<p style='color:#888;font-style:italic;'>No statistics for this module.</p>"
    cols = ["metric", "n_pairs", "median_primary", "median_recurrent",
            "median_delta", "mean_delta", "wilcoxon_stat", "p_value",
            "p_value_adj_bh", "note"]
    th = "".join(f"<th>{c}</th>" for c in cols)

    def fmt(v: object) -> str:
        if isinstance(v, float):
            if pd.isna(v):
                return "—"
            return f"{v:.4g}"
        return str(v)

    body = "".join(
        "<tr>" + "".join(f"<td>{fmt(r[c])}</td>" for c in cols) + "</tr>"
        for _, r in sub.iterrows())
    return f"<table><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table>"


def _df_table_html(df: pd.DataFrame, empty: str = "No data.") -> str:
    """Render a full DataFrame as an HTML table (4 sig-fig floats)."""
    if df is None or df.empty:
        return f"<p style='color:#888;font-style:italic;'>{empty}</p>"
    th = "".join(f"<th>{c}</th>" for c in df.columns)

    def fmt(v: object) -> str:
        if isinstance(v, float):
            if pd.isna(v):
                return "—"
            return f"{v:.4g}"
        return str(v)

    body = "".join(
        "<tr>" + "".join(f"<td>{fmt(r[c])}</td>" for c in df.columns) + "</tr>"
        for _, r in df.iterrows())
    return f'<div style="overflow-x:auto;">' \
           f'<table><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table></div>'


def generate_report(cohort: pd.DataFrame, fs_cohort: pd.DataFrame,
                    summary: pd.DataFrame, paired: pd.DataFrame,
                    landscape: pd.DataFrame, out_dir: Path, n_input_samples: int,
                    mod3: pd.DataFrame, mod4: pd.DataFrame, mod5: pd.DataFrame,
                    stage3_stats: pd.DataFrame, module_imgs: Dict[str, str],
                    overlap_data: Dict[str, dict]):
    out_dir.mkdir(parents=True, exist_ok=True)

    if cohort.empty:
        (out_dir / "cohort_report.html").write_text(
            f"<!DOCTYPE html><html><head><style>{CSS}</style></head>"
            "<body><div class='page'><h1>NMD Cohort Report</h1>"
            "<p class='sub'>No candidates found across cohort.</p>"
            "</div></body></html>")
        print("[REPORT] Cohort empty — wrote stub")
        return

    # Static plots
    img_landscape = plot_variant_landscape(cohort)
    img_tiers     = plot_tiers(fs_cohort) if not fs_cohort.empty else ""
    top_for_ic50  = (fs_cohort.sort_values("best_ic50").head(50)
                     if "best_ic50" in fs_cohort.columns else fs_cohort.head(0))
    img_ic50      = plot_ic50(top_for_ic50)
    img_method    = plot_nmd_breakdown(fs_cohort) if not fs_cohort.empty else ""
    img_conf      = plot_confidence(fs_cohort) if not fs_cohort.empty else ""
    img_t_v_m     = plot_tier_by_timepoint(fs_cohort)
    img_paired    = plot_paired(paired)
    img_genes     = plot_top_genes(fs_cohort)

    hla_df = hla_allele_breakdown(fs_cohort) if not fs_cohort.empty else pd.DataFrame()

    # Top tables (FS only)
    t1 = fs_cohort[fs_cohort["priority_tier"] == "TIER1"]
    t2 = fs_cohort[fs_cohort["priority_tier"] == "TIER2"]
    t3 = fs_cohort[fs_cohort["priority_tier"] == "TIER3_control"]
    if "best_ic50" in fs_cohort.columns:
        t1 = t1.sort_values("best_ic50").head(30)
        t2 = t2.sort_values("best_ic50").head(30)
        t3 = t3.sort_values("best_ic50").head(30)

    COLS = ["sample","patient","timepoint","Gene Name","hla_allele",
            "MT Epitope Seq","best_ic50","best_ic50_method","nmd_consensus",
            "nmd_confidence","nmd_confidence_score","nmd_rule_explanation"]
    COLS = [c for c in COLS if c in fs_cohort.columns]

    # Cards
    n_smp_loaded = cohort["sample"].nunique()
    n_pat        = cohort["patient"].nunique()
    n_total      = len(cohort)
    n_fs         = len(fs_cohort)
    n_t1         = int((fs_cohort["priority_tier"]=="TIER1").sum()) if not fs_cohort.empty else 0
    n_t2         = int((fs_cohort["priority_tier"]=="TIER2").sum()) if not fs_cohort.empty else 0
    n_s          = int((fs_cohort["nmd_consensus"]=="SENSITIVE").sum()) if not fs_cohort.empty else 0
    n_i          = int((fs_cohort["nmd_consensus"]=="INSENSITIVE").sum()) if not fs_cohort.empty else 0
    n_dis        = int((fs_cohort["nmd_confidence"]=="methods_disagree").sum()) if not fs_cohort.empty else 0
    n_no_data    = int((fs_cohort["nmd_confidence"]=="no_data").sum()) if not fs_cohort.empty else 0
    n_unclass    = int((fs_cohort["priority_tier"]=="UNCLASSIFIED").sum()) if not fs_cohort.empty else 0

    n_total_pat = 28  # full cohort patient count from sample list (constant for this study)

    # Per-patient T+M coverage (answers "why 27 patients but 49 samples":
    # both timepoints have candidates for some patients, only one for others,
    # neither for patient 5 who has very low purity in both samples).
    if cohort.empty:
        n_both = n_one = n_none = 0
    else:
        cov = cohort.groupby("patient")["timepoint"].nunique()
        n_both = int((cov == 2).sum())
        n_one  = int((cov == 1).sum())
        n_none = int(n_total_pat - n_both - n_one)

    def conf_cards():
        return ('<div style="display:flex;gap:12px;flex-wrap:wrap;margin-bottom:16px;">' +
                "".join(f'<div class="card"><div class="cl">Score {s}</div>'
                        f'<div class="cv">{int((fs_cohort["nmd_confidence_score"]==s).sum()) if not fs_cohort.empty else 0}</div></div>'
                        for s in [3,2,1,0]) + "</div>")

    def hla_tbl():
        if hla_df is None or hla_df.empty:
            return "<p style='color:#888;font-style:italic;'>No Tier 1 or Tier 2 candidates.</p>"
        cols = ["HLA Allele","Candidates","Tier 1","Best IC50 (nM)","Median IC50 (nM)"]
        th = "".join(f"<th>{c}</th>" for c in cols)
        body = "".join(
            f"<tr><td>{r.hla_allele}</td><td>{r.n_candidates}</td><td>{r.n_tier1}</td>"
            f"<td>{r.best_ic50:.1f}</td><td>{r.median_ic50:.1f}</td></tr>"
            for _, r in hla_df.iterrows())
        return f"<table><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table>"

    patient_data_json = json.dumps(build_patient_data(fs_cohort), separators=(",",":"))
    widget_block = WIDGET_HTML.replace("__PATIENT_DATA__", patient_data_json)

    overlap_json = json.dumps(overlap_data, separators=(",", ":"))
    overlap_block = OVERLAP_WIDGET_HTML.replace("__OVERLAP_DATA__", overlap_json)

    # ── Module 3/4/5 sections ────────────────────────────────────────────────
    module3_block = f"""
<h2>Module 3 — NMD-sensitive mutation enrichment</h2>
<p class="sub">
  Variant-level NMD aggregation over truncating (frameshift) variants. Each
  variant is classified once (worst-case across its peptides). NMD fraction =
  sensitive / (sensitive + insensitive). Counts are deduplicated to the variant
  level.
</p>
{_img(module_imgs.get("m3_paired",""), "M3 Fig 1. Paired NMD-sensitive fraction (FS), primary vs recurrent.")}
{_img(module_imgs.get("m3_delta",""), "M3 Fig 2. ΔNMD fraction (recurrent − primary), sorted.")}
{_img(module_imgs.get("m3_stacked",""), "M3 Fig 3. NMD-sensitive vs escape FS variants per sample.")}
{_img(module_imgs.get("m3_scatter",""), "M3 Fig 4. Truncating burden vs NMD-sensitive count, coloured by timepoint.")}
<h3>Per-sample NMD summary</h3>
{_df_table_html(mod3, "No Module 3 data.")}
<h3>Module 3 paired statistics</h3>
{_stats_table_html(stage3_stats, "3")}
"""

    module4_block = f"""
<h2>Module 4 — Neoepitope burden</h2>
<p class="sub">
  Per-sample binder burden across all candidate peptides. Strong &lt;50 nM,
  Weak 50–500 nM, Non-binder ≥500 nM. Clonality from Tumor DNA VAF
  (Clonal ≥0.30). fraction_clonal = clonal / (clonal + subclonal) over binders.
</p>
{_img(module_imgs.get("m4_paired",""), "M4 Fig 1. Paired total and strong binder burden.")}
{_img(module_imgs.get("m4_stacked",""), "M4 Fig 2. Clonal vs subclonal binders per sample.")}
{_img(module_imgs.get("m4_scatter",""), "M4 Fig 3. TMB (Stage 1, log) vs strong-binder burden.") if module_imgs.get("m4_scatter") else "<p class='note'>TMB-vs-neoepitope scatter skipped — Stage 1 summary not provided.</p>"}
{_img(module_imgs.get("m4_delta",""), "M4 Fig 4. Δ strong binders (recurrent − primary), sorted.")}
<h3>Per-sample neoepitope summary</h3>
{_df_table_html(mod4, "No Module 4 data.")}
<h3>Module 4 paired statistics</h3>
{_stats_table_html(stage3_stats, "4")}
"""

    module5_block = f"""
<h2>Module 5 — Neoepitopes by NMD <span style="font-size:12px;color:#888;">(the central mechanistic result)</span></h2>
<p class="sub">
  Neoepitopes (Strong/Weak binders) split by their variant's NMD class.
  fraction_escape = escape / (escape + sensitive); UNCERTAIN/UNKNOWN excluded.
  The headline figure tests whether recurrence shifts neoepitopes toward NMD
  escape (immune visibility) or toward NMD-sensitive silencing.
</p>
<div style="border:2px solid #D85A30;border-radius:8px;padding:8px 12px;margin:12px 0;background:#fffaf7;">
  <h3 style="margin-top:4px;color:#9a3a14;">Headline figure — paired NMD-escape fraction (strong binders)</h3>
  {_img(module_imgs.get("m5_paired",""), "M5 Fig 1 (HEADLINE). Paired fraction_escape_strong, primary vs recurrent, with p-value and median Δ.")}
</div>
{_img(module_imgs.get("m5_stacked",""), "M5 Fig 2. NMD-sensitive vs escape strong neoepitopes per sample.")}
{_img(module_imgs.get("m5_delta",""), "M5 Fig 3. Δ fraction_escape_strong (recurrent − primary), sorted.")}
{_img(module_imgs.get("m5_paired_count",""), "M5 Fig 4. Paired NMD-sensitive neoepitope count (strong).")}
<h3>Per-sample neoepitope-by-NMD summary</h3>
{_df_table_html(mod5, "No Module 5 data.")}
<h3>Module 5 paired statistics</h3>
{_stats_table_html(stage3_stats, "5")}
"""

    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
<title>GBM NMD Cohort Report</title><style>{CSS}</style></head>
<body><div class="page">

<h1>GBM Pipeline — Stage 3: NMD Cohort Summary</h1>
<p style="color:#888;font-size:13px;margin-top:4px;">
  GBM NMD-Neoantigen Pipeline &middot; github.com/paleslui/gbm-nmd-pipeline
</p>

<div class="intro-box">
  This report summarises NMD-actionable neoantigen candidates across the cohort.
  The hypothesis is that TMZ-induced frameshift mutations create premature
  termination codons (PTCs) that are silenced by nonsense-mediated mRNA decay
  (NMD); inhibiting NMD could expose these neoantigens to immune recognition.
  Frameshift variants are the only neoantigen class that creates novel
  peptides actionable by NMD inhibition; missense variants are not NMD-relevant
  and stop-gained variants do not produce peptide neoantigens. The variant
  landscape below is shown for context, but the rest of this report focuses
  exclusively on frameshift-derived candidates.
</div>

{METHODS_BLOCK}

<h2>Cohort overview</h2>
{_card("Samples", f"{n_smp_loaded}/{n_input_samples}", "≥1 candidate")}
{_card("Patients", f"{n_pat}/{n_total_pat}", "≥1 candidate")}
{_card("Patient T+M", f"{n_both}/{n_one}/{n_none}", "both / one only / neither")}
{_card("Total candidates", n_total, "all variant types")}
{_card("FS candidates", n_fs, "NMD-actionable")}
{_card("TIER 1", n_t1, "FS, IC50&lt;50 nM")}
{_card("TIER 2", n_t2)}
{_card("NMD-sensitive", n_s, "FS only")}
{_card("NMD-insensitive", n_i, "FS only")}
{_card("Unclassified", n_unclass, f"{n_no_data} no data · {n_dis} methods disagree")}
<p class="sub" style="margin-top:14px;">
  Note on samples vs patients: the cohort contains paired primary (T) and recurrent (M) samples,
  so {n_total_pat} patients = {n_total_pat*2} samples. Of the {n_input_samples - n_smp_loaded} samples with 0 candidates,
  most are paired with a non-empty timepoint in the same patient (so the patient still counts as represented).
  Only when both T and M timepoints are empty does a patient drop out of the cohort entirely.
</p>

<h2>Variant landscape — input candidates by type</h2>
<p class="sub">
  Distribution of all <strong>{n_total}</strong> pVACseq candidates across the cohort.
  Frameshift (FS) is the only NMD-actionable variant type for neoantigen-based
  immunotherapy. The remaining sections of this report show NMD analysis on the
  <strong>{n_fs} FS candidates</strong> only.
</p>
{_img(img_landscape, f"Fig 1. Variant type distribution. Of {n_total} candidates, {n_fs} are frameshift (NMD-actionable).")}
{_tbl(landscape, list(landscape.columns), "No variant data.")}
<div class="note">
  <strong>Why we focus on frameshift only:</strong>
  <em>Missense</em> variants are not NMD-relevant (they replace one amino acid without
  introducing a PTC). <em>Stop-gained</em> variants are NMD-relevant but truncate
  the protein without generating novel peptide sequences — pVACseq does not
  emit peptide candidates for them, so they are not actionable for neoantigen
  vaccines. <em>Inframe</em> indels do not introduce PTCs.
</div>

<h2>Cohort summary table</h2>
{_tbl(summary, list(summary.columns), "Empty cohort")}

<h2>Tier distribution (FS only)</h2>
{_img(img_tiers, "Fig 2. Cohort-wide candidate counts by priority tier (frameshift only).")}

<h2>Tier distribution: primary (T) vs recurrent (M)</h2>
{_img(img_t_v_m, "Fig 3. Tier counts split by timepoint. Recurrent enrichment of TIER1 candidates is the core thesis hypothesis.")}

<h2>Per-patient paired comparison</h2>
{_img(img_paired, "Fig 4. TIER1 candidate count per patient, primary vs recurrent. Sorted by ΔTIER1 (recurrent − primary), patients gaining the most from T→M first.")}
<p class="sub">Per-patient counts (sorted by ΔTIER1, then ΔTotal):</p>
{_tbl(paired, list(paired.columns), "No paired data.")}

{module3_block}

{module4_block}

{module5_block}

{overlap_block}

{widget_block}

<h2>IC50 distribution by NMD class — top 50 strongest binders (FS only)</h2>
{_img(img_ic50, "Fig 5. Top 50 frameshift candidates across the cohort, sorted by best IC50, colored by NMD consensus.")}

<h2>Top genes producing TIER1 candidates</h2>
{_img(img_genes, "Fig 6. Genes contributing the most TIER1 candidates across the cohort (frameshift only).")}

<h2>Per-HLA allele breakdown (Tier 1 + Tier 2)</h2>
<p class="sub">Tier 1 and Tier 2 candidates across HLA alleles. Alleles with multiple
high-confidence binders are the strongest therapeutic targets.</p>
{hla_tbl()}

<h2>NMD classification per method (cohort-wide, FS only)</h2>
{_img(img_method, "Fig 7. Cohort-wide NMD classification per scoring method.")}

<h2>Confidence score distribution (cohort-wide, FS only)</h2>
{_img(img_conf, "Fig 8. Cohort-wide ensemble confidence distribution.")}
<p class="sub">Score 3 = both methods agree; 2 = single method; 1 = disagree; 0 = no data.</p>
{conf_cards()}

<h2>TIER 1 — NMD-sensitive + IC50 &lt; 50 nM (top 30 across cohort by IC50)</h2>
{_tbl(t1, COLS, "No Tier 1 candidates in this cohort.")}

<h2>TIER 2 — NMD-sensitive + IC50 50–500 nM (top 30 across cohort by IC50)</h2>
{_tbl(t2, COLS, "No Tier 2 candidates in this cohort.")}

<h2>TIER 3 — NMD-insensitive controls (top 30 across cohort by IC50)</h2>
{_tbl(t3, COLS, "No Tier 3 candidates in this cohort.")}

</div></body></html>"""

    out = out_dir / "cohort_report.html"
    out.write_text(html, encoding="utf-8")
    print(f"[REPORT] Saved {out}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description="GBM Stage 3 — NMD cohort summary")
    ap.add_argument("--input_dir", required=True, type=Path,
                    help="Directory containing per-sample subdirs (each with nmd_scored_candidates.tsv)")
    ap.add_argument("--out_dir", required=True, type=Path,
                    help="Where cohort_*.{tsv,html} are written")
    ap.add_argument("--stage1_summary", type=Path, default=None,
                    help="Stage 1 sample_mutation_summary.tsv (for truncating_total / TMB joins)")
    ap.add_argument("--paired_overlap", type=Path, default=None,
                    help="Stage 1 paired_variant_overlap.tsv (for the overlap widget)")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}\n  GBM Pipeline — Stage 3 cohort summary\n{'='*60}\n")
    print(f"[INFO] Reading per-sample dirs from: {args.input_dir}")

    # ── optional Stage 1 joins (graceful fallback) ──────────────────────────
    stage1: Optional[pd.DataFrame] = None
    if args.stage1_summary and args.stage1_summary.is_file():
        stage1 = pd.read_csv(args.stage1_summary, sep="\t")
        stage1["sample"] = stage1["sample"].astype(str)
        print(f"[INFO] Stage 1 summary loaded: {len(stage1)} samples")
    elif args.stage1_summary:
        print(f"[WARN] Stage 1 summary not found at {args.stage1_summary}; "
              "stage1-dependent columns will be NaN")
    else:
        print("[INFO] No --stage1_summary given; stage1-dependent columns will be NaN")

    overlap_df: Optional[pd.DataFrame] = None
    if args.paired_overlap and args.paired_overlap.is_file():
        overlap_df = pd.read_csv(args.paired_overlap, sep="\t")
        print(f"[INFO] Paired overlap loaded: {len(overlap_df)} patients")
    elif args.paired_overlap:
        print(f"[WARN] Paired overlap not found at {args.paired_overlap}; "
              "overlap widget will be empty")

    n_input_samples = sum(1 for p in args.input_dir.iterdir() if p.is_dir() and SAMPLE_RE.match(p.name))
    cohort = load_cohort(args.input_dir)

    # Augment candidates with binder_class + clonality (appended columns).
    cohort = augment_candidates(cohort)

    fs_cohort = cohort[cohort["Variant Type"] == "FS"].copy() if not cohort.empty else cohort
    print(f"[INFO] Filtered to FS-only: {len(fs_cohort)} candidates from {fs_cohort['sample'].nunique() if not fs_cohort.empty else 0} samples")

    landscape = variant_landscape(cohort)
    summary   = cohort_summary(cohort, fs_cohort, n_input_samples)
    paired    = per_patient_paired(fs_cohort)

    # ── Max Modules 3/4/5 ────────────────────────────────────────────────────
    variant_tbl = variant_nmd_table(cohort)
    mod3 = module3_nmd_summary(cohort, variant_tbl, stage1) if not cohort.empty else pd.DataFrame()
    mod4 = module4_neoepitope_summary(cohort) if not cohort.empty else pd.DataFrame()
    mod5 = module5_neoepitope_nmd_summary(cohort, variant_tbl) if not cohort.empty else pd.DataFrame()

    # Wide per-sample frame for paired stats + plots.
    have_stage1 = stage1 is not None
    if not cohort.empty:
        wide = mod3.merge(mod4, on=["sample", "patient", "timepoint"]) \
                   .merge(mod5, on=["sample", "patient", "timepoint"])
        if have_stage1:
            tmb = stage1[["sample"]].copy()
            for col in ("tmb_per_mb", "sbs11_pct", "snv", "indel", "coding"):
                if col in stage1.columns:
                    tmb[col] = stage1[col]
            wide = wide.merge(tmb, on="sample", how="left")
        stage3_stats = compute_paired_stats_stage3(wide)
        module_imgs = make_module_plots(wide, stage3_stats, args.out_dir, have_stage1)
    else:
        wide = pd.DataFrame()
        stage3_stats = pd.DataFrame(columns=["module", "metric", "n_pairs",
                                             "median_primary", "median_recurrent",
                                             "median_delta", "mean_delta",
                                             "wilcoxon_stat", "p_value",
                                             "p_value_adj_bh", "note"])
        module_imgs = {}

    overlap_data = build_overlap_data(overlap_df, variant_tbl, stage1)

    # ── write TSVs (existing + new; never delete existing outputs) ───────────
    if not cohort.empty:
        cohort.to_csv(args.out_dir / "cohort_candidates.tsv", sep="\t", index=False)
        if not fs_cohort.empty:
            fs_cohort[fs_cohort["priority_tier"]=="TIER1"].to_csv(
                args.out_dir / "cohort_tier1.tsv", sep="\t", index=False)
    summary.to_csv  (args.out_dir / "cohort_summary.tsv",  sep="\t", index=False)
    paired.to_csv   (args.out_dir / "cohort_paired.tsv",   sep="\t", index=False)
    landscape.to_csv(args.out_dir / "cohort_landscape.tsv", sep="\t", index=False)
    mod3.to_csv(args.out_dir / "per_sample_nmd_summary.tsv", sep="\t", index=False)
    mod4.to_csv(args.out_dir / "per_sample_neoepitope_summary.tsv", sep="\t", index=False)
    mod5.to_csv(args.out_dir / "per_sample_neoepitope_nmd_summary.tsv", sep="\t", index=False)
    stage3_stats.to_csv(args.out_dir / "paired_stats_stage3.tsv", sep="\t", index=False)

    generate_report(cohort, fs_cohort, summary, paired, landscape, args.out_dir,
                    n_input_samples, mod3, mod4, mod5, stage3_stats, module_imgs,
                    overlap_data)

    print(f"\n{'='*60}\n  Cohort summary done. Output: {args.out_dir.resolve()}\n{'='*60}\n")


if __name__ == "__main__":
    main()
