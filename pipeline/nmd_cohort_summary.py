#!/usr/bin/env python3
"""
nmd_cohort_summary.py — GBM Pipeline Stage 3: NMD Cohort Summary
Part of the GBM NMD-Neoantigen Pipeline.
https://github.com/paleslui/gbm-nmd-pipeline

Aggregates per-sample NMD scoring outputs (Stage 3 per_sample/) into:
  - cohort_candidates.tsv : all candidates with sample/patient/timepoint
  - cohort_summary.tsv    : cohort-level counts (tier, NMD, by patient)
  - cohort_tier1.tsv      : TIER1 high-priority candidates only
  - cohort_paired.tsv     : per-patient T (primary) vs M (recurrent) comparison
  - cohort_report.html    : single HTML report mirroring the per-sample style
                            (same plotting functions imported from nmd_scoring.py)
                            plus cohort-specific plots:
                              - tier-by-timepoint
                              - per-patient paired bar chart
                              - top genes producing TIER1 candidates

Usage:
  python nmd_cohort_summary.py --input_dir <per_sample_dir> --out_dir <cohort_dir>

  --input_dir : dir containing one subdir per sample (each with
                nmd_scored_candidates.tsv); typically run_*/3_nmd_analysis/per_sample/
  --out_dir   : where cohort_*.{tsv,html} are written
"""
import argparse, re, sys
from pathlib import Path

import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Reuse plotting + style from per-sample scorer for consistent look
sys.path.insert(0, str(Path(__file__).parent))
from nmd_scoring import (
    plot_tiers, plot_ic50, plot_nmd_breakdown, plot_confidence,
    hla_allele_breakdown, _b64fig, COL_S, COL_I, COL_U,
)

# Sample dir name like "11_T" or "11_M"
SAMPLE_RE = re.compile(r"^(\d+)_([TM])$")


# ═════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═════════════════════════════════════════════════════════════════════════════

def load_cohort(per_sample_dir: Path) -> pd.DataFrame:
    """Load all per-sample TSVs from per_sample_dir/{sample}/nmd_scored_candidates.tsv"""
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
    # Coerce numeric cols used downstream (in nmd_scoring.py these are floats)
    for col in ("best_ic50", "median_ic50", "nmd_confidence_score"):
        if col in cohort.columns:
            cohort[col] = pd.to_numeric(cohort[col], errors="coerce")
    print(f"[INFO] Cohort total: {len(cohort)} candidates from {cohort['sample'].nunique()} samples")
    return cohort


# ═════════════════════════════════════════════════════════════════════════════
# AGGREGATE TABLES
# ═════════════════════════════════════════════════════════════════════════════

def cohort_summary(cohort: pd.DataFrame, n_input_samples: int) -> pd.DataFrame:
    """Cohort-level aggregate counts."""
    if cohort.empty:
        return pd.DataFrame([("No candidates", 0)], columns=["Metric", "Count"])
    rows = [
        ("Total candidates",                  len(cohort)),
        ("Samples with ≥1 candidate",         cohort["sample"].nunique()),
        ("Samples with 0 candidates",         n_input_samples - cohort["sample"].nunique()),
        ("Patients represented",              cohort["patient"].nunique()),
        ("", ""),
        ("NMD-SENSITIVE",                     int((cohort["nmd_consensus"]=="SENSITIVE").sum())),
        ("NMD-INSENSITIVE",                   int((cohort["nmd_consensus"]=="INSENSITIVE").sum())),
        ("UNCERTAIN / UNKNOWN",               int(cohort["nmd_consensus"].isin(["UNCERTAIN","UNKNOWN"]).sum())),
        ("NOT_APPLICABLE",                    int((cohort.get("nmd_rules", pd.Series())=="NOT_APPLICABLE").sum())),
        ("", ""),
        ("TIER1 (NMD-sensitive, IC50<50)",    int((cohort["priority_tier"]=="TIER1").sum())),
        ("TIER2 (NMD-sensitive, IC50<500)",   int((cohort["priority_tier"]=="TIER2").sum())),
        ("TIER3 controls (NMD-insensitive)",  int((cohort["priority_tier"]=="TIER3_control").sum())),
        ("Unclassified",                      int((cohort["priority_tier"]=="UNCLASSIFIED").sum())),
        ("", ""),
        ("High confidence (3 — both agree)",  int((cohort["nmd_confidence_score"]==3).sum())),
        ("Medium confidence (2 — single)",    int((cohort["nmd_confidence_score"]==2).sum())),
        ("Low confidence (1 — disagree)",     int((cohort["nmd_confidence_score"]==1).sum())),
        ("No data (0)",                       int((cohort["nmd_confidence_score"]==0).sum())),
    ]
    return pd.DataFrame(rows, columns=["Metric", "Count"])


def per_patient_paired(cohort: pd.DataFrame) -> pd.DataFrame:
    """T vs M comparison per patient: counts of candidates and tiers."""
    if cohort.empty:
        return pd.DataFrame()
    rows = []
    for patient, sub in cohort.groupby("patient"):
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
# (per-sample plots are imported from nmd_scoring.py for visual consistency)
# ═════════════════════════════════════════════════════════════════════════════

COL_T = "#378ADD"  # primary
COL_M = "#D85A30"  # recurrent

def plot_tier_by_timepoint(cohort: pd.DataFrame) -> str:
    """Stacked-grouped: TIER1/2/3 distribution split primary vs recurrent."""
    if cohort.empty:
        return ""
    counts = (cohort.groupby(["timepoint","priority_tier"]).size()
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
    ax.set_ylabel("Candidates")
    ax.set_title("Tier distribution: primary (T) vs recurrent (M)", fontsize=12)
    ax.legend(fontsize=9)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    return _b64fig(fig)


def plot_paired(paired: pd.DataFrame) -> str:
    """Per-patient TIER1 count: T vs M side-by-side bars."""
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
    ax.set_title("TIER1 NMD-sensitive neoantigens per patient: primary vs recurrent",
                 fontsize=12)
    ax.legend(fontsize=9)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    return _b64fig(fig)


def plot_top_genes(cohort: pd.DataFrame, n: int = 20) -> str:
    """Top genes producing TIER1 candidates (horizontal bar)."""
    if cohort.empty:
        return ""
    tier1 = cohort[cohort["priority_tier"] == "TIER1"]
    if tier1.empty:
        return ""
    counts = tier1["Gene Name"].value_counts().head(n)
    fig, ax = plt.subplots(figsize=(8, max(3, 0.32 * len(counts))))
    counts.iloc[::-1].plot(kind="barh", ax=ax, color="#F0C040")
    for i, v in enumerate(counts.iloc[::-1]):
        ax.text(v + 0.05, i, str(int(v)), va="center", fontsize=9)
    ax.set_xlabel("TIER1 candidate count")
    ax.set_title(f"Top {len(counts)} genes producing TIER1 candidates", fontsize=12)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    return _b64fig(fig)


# ═════════════════════════════════════════════════════════════════════════════
# HTML REPORT
# Mirrors the per-sample report style (CSS, sections, plot calls).
# Inlines explanatory text equivalent to nmd_scoring.py's "What is this report?",
# "NMD Scoring Methods", "Priority Tiers" blocks so a reader of the cohort
# report alone has full context.
# ═════════════════════════════════════════════════════════════════════════════

CSS = ("*{box-sizing:border-box} body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',"
       "sans-serif;font-size:14px;line-height:1.6;color:#1a1a1a;background:#fff;margin:0;padding:0}"
       ".page{max-width:1200px;margin:0 auto;padding:40px 32px}"
       "h1{font-size:24px;font-weight:500;margin-bottom:4px}"
       "h2{font-size:17px;font-weight:500;border-bottom:1px solid #e0e0e0;padding-bottom:6px;margin:36px 0 8px}"
       "h3{font-size:14px;font-weight:500;margin:16px 0 6px;color:#444}"
       "p.sub{color:#666;font-size:13px;margin:0 0 16px;line-height:1.5}"
       "table{border-collapse:collapse;width:100%;font-size:12px;margin-top:8px}"
       "th{background:#f4f4f4;padding:6px 10px;text-align:left;font-weight:500;border-bottom:2px solid #ddd}"
       "td{padding:6px 10px;border-bottom:1px solid #eee}"
       ".card{background:#f8f8f8;border-radius:8px;padding:14px 18px;display:inline-block;min-width:140px;margin:4px}"
       ".cv{font-size:22px;font-weight:500} .cl{font-size:11px;color:#777}"
       ".note{background:#fffbe6;border-left:3px solid #f0c040;padding:10px 14px;margin:12px 0;font-size:13px;color:#555}")


def _card(label, val, color=""):
    style = "" if not color else f' style="color:{color}"'
    return (f'<div class="card"><div class="cl">{label}</div>'
            f'<div class="cv"{style}>{val}</div></div>')


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


def generate_report(cohort: pd.DataFrame, summary: pd.DataFrame, paired: pd.DataFrame,
                    out_dir: Path, n_input_samples: int):
    out_dir.mkdir(parents=True, exist_ok=True)

    # If cohort is empty, write a stub
    if cohort.empty:
        (out_dir / "cohort_report.html").write_text(
            f"<!DOCTYPE html><html><head><style>{CSS}</style></head>"
            "<body><div class='page'><h1>NMD Cohort Report</h1>"
            "<p class='sub'>No candidates found across cohort.</p>"
            "</div></body></html>")
        print("[REPORT] Cohort empty — wrote stub")
        return

    # Generate plots (per-sample-style, applied to cohort)
    img_tiers   = plot_tiers(cohort)
    # IC50 scatter is one row per candidate; cap at top-50 to stay readable
    top_for_ic50 = (cohort.sort_values("best_ic50").head(50)
                    if "best_ic50" in cohort.columns else cohort.head(0))
    img_ic50    = plot_ic50(top_for_ic50)
    img_method  = plot_nmd_breakdown(cohort)
    img_conf    = plot_confidence(cohort)
    # Cohort-specific
    img_t_v_m   = plot_tier_by_timepoint(cohort)
    img_paired  = plot_paired(paired)
    img_genes   = plot_top_genes(cohort)

    hla_df = hla_allele_breakdown(cohort)

    # Subsets for tables
    t1 = cohort[cohort["priority_tier"] == "TIER1"]
    t2 = cohort[cohort["priority_tier"] == "TIER2"]
    t3 = cohort[cohort["priority_tier"] == "TIER3_control"]

    # Show top 30 of each tier sorted by best_ic50 (ascending = strongest binders first)
    if "best_ic50" in cohort.columns:
        t1 = t1.sort_values("best_ic50").head(30)
        t2 = t2.sort_values("best_ic50").head(30)
        t3 = t3.sort_values("best_ic50").head(30)

    COLS = ["sample","patient","timepoint","Gene Name","Variant Type","hla_allele",
            "MT Epitope Seq","best_ic50","best_ic50_method","nmd_consensus",
            "nmd_confidence","nmd_confidence_score","nmd_rule_explanation"]
    COLS = [c for c in COLS if c in cohort.columns]

    # Counters for headline cards
    n_tot = len(cohort)
    n_s   = int((cohort["nmd_consensus"]=="SENSITIVE").sum())
    n_i   = int((cohort["nmd_consensus"]=="INSENSITIVE").sum())
    n_u   = int(cohort["nmd_consensus"].isin(["UNCERTAIN","UNKNOWN"]).sum())
    n_t1  = int((cohort["priority_tier"]=="TIER1").sum())
    n_t2  = int((cohort["priority_tier"]=="TIER2").sum())
    n_dis = int((cohort["nmd_confidence"]=="methods_disagree").sum())
    n_pat = cohort["patient"].nunique()
    n_smp = cohort["sample"].nunique()

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

    def conf_cards():
        return ('<div style="display:flex;gap:12px;flex-wrap:wrap;margin-bottom:16px;">' +
                "".join(f'<div class="card"><div class="cl">Score {s}</div>'
                        f'<div class="cv">{int((cohort["nmd_confidence_score"]==s).sum())}</div></div>'
                        for s in [3,2,1,0]) + "</div>")

    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
<title>GBM NMD Cohort Report</title><style>{CSS}</style></head>
<body><div class="page">

<h1>GBM Pipeline — Stage 3: NMD Cohort Summary</h1>
<p style="color:#888;font-size:13px;margin-top:4px;">
  GBM NMD-Neoantigen Pipeline &middot; github.com/paleslui/gbm-nmd-pipeline
</p>

<h2>What is this report?</h2>
<p class="sub">
  Aggregates per-sample NMD scoring results across the whole cohort.
  Each pVACseq neoantigen candidate has been classified by NMD sensitivity using an
  ensemble of two methods (VEP NMD plugin + Lindeboom et al. 2019 rules).
  The core hypothesis is that TMZ-induced frameshift mutations create PTCs silenced by NMD —
  NMD inhibition could expose these neoantigens to immune recognition. This report
  compares primary (T) and recurrent (M) tumours per patient and identifies the
  highest-priority candidates across the cohort.
</p>

<h2>Cohort overview</h2>
{_card("Samples (loaded)", n_smp)}
{_card("Samples (input)", n_input_samples)}
{_card("Patients", n_pat)}
{_card("Total candidates", n_tot)}
{_card("NMD-sensitive", n_s, COL_S)}
{_card("NMD-insensitive", n_i, COL_I)}
{_card("Uncertain / N/A", n_u, COL_U)}
{_card("Tier 1", n_t1)}
{_card("Tier 2", n_t2)}
{_card("Methods disagree", n_dis)}

<h2>Cohort summary table</h2>
{_tbl(summary, list(summary.columns), "Empty cohort")}

<h2>Tier distribution (cohort-wide)</h2>
{_img(img_tiers, "Fig 1. Cohort-wide candidate counts by priority tier.")}

<h2>Tier distribution: primary (T) vs recurrent (M)</h2>
{_img(img_t_v_m, "Fig 2. Tier counts split by timepoint. Recurrent enrichment of TIER1 candidates is the core thesis hypothesis.")}

<h2>Per-patient paired comparison</h2>
{_img(img_paired, "Fig 3. TIER1 candidate count per patient, primary vs recurrent. Sorted by ΔTIER1 (recurrent - primary), patients gaining the most from T→M first.")}
<p class="sub">Per-patient counts (sorted by ΔTIER1, then ΔTotal):</p>
{_tbl(paired, list(paired.columns), "No paired data.")}

<h2>IC50 distribution by NMD class — top 50 strongest binders</h2>
{_img(img_ic50, "Fig 4. Top 50 candidates across the cohort, sorted by best IC50, colored by NMD consensus.")}

<h2>Top genes producing TIER1 candidates</h2>
{_img(img_genes, "Fig 5. Genes contributing the most TIER1 candidates across the cohort.")}

<h2>Per-HLA allele breakdown</h2>
<p class="sub">Tier 1 and Tier 2 candidates across HLA alleles. Alleles with multiple
high-confidence binders are the strongest therapeutic targets.</p>
{hla_tbl()}

<h2>NMD classification per method (cohort-wide)</h2>
{_img(img_method, "Fig 6. Cohort-wide NMD classification per scoring method.")}

<h2>Confidence score distribution (cohort-wide)</h2>
{_img(img_conf, "Fig 7. Cohort-wide ensemble confidence distribution.")}
<p class="sub">Score 3 = both methods agree; 2 = single method; 1 = disagree; 0 = no data.</p>
{conf_cards()}

<h2>NMD scoring methods</h2>
<h3>Method 1 — VEP NMD plugin</h3>
<p class="sub">
  The VEP NMD plugin (Ensembl v105+) annotates truncating variants with
  <code>NMD_escaping_variant</code> when the PTC is predicted to escape NMD. An empty
  NMD field on a truncating variant means NMD is triggered (SENSITIVE).
</p>

<h3>Method 2 — Lindeboom rules (Nat Cell Biol 2019)</h3>
<p class="sub">Applied to frameshift and stop-gained variants. Rules in priority order:</p>
<ul style="font-size:13px;color:#444;margin:0 0 12px 20px;line-height:1.8">
  <li><strong>Rule 4 — Start-proximal PTC (&lt;150nt):</strong> Pioneer round completes before NMD surveillance → INSENSITIVE</li>
  <li><strong>Rule 1 — Last exon:</strong> No downstream EJC to trigger NMD → INSENSITIVE</li>
  <li><strong>Rule 3 — Long exon (&gt;407nt):</strong> EJC too far downstream → INSENSITIVE</li>
  <li><strong>Rule 2 — 55nt boundary:</strong> PTC &gt;55nt upstream of last EJC → SENSITIVE (canonical NMD)</li>
</ul>
<div class="note">⚠ NMD rules only apply to truncating mutations. Missense variants
(NOT_APPLICABLE) are not subject to NMD but remain immunogenic neoantigen candidates.</div>

<h2>Priority tiers</h2>
<ul style="font-size:13px;color:#444;margin:0 0 12px 20px;line-height:1.8">
  <li><strong>Tier 1 — NMD-sensitive + IC50 &lt;50 nM:</strong> Primary therapeutic targets — NMD inhibition could expose these neoantigens.</li>
  <li><strong>Tier 2 — NMD-sensitive + IC50 50–500 nM:</strong> Moderate binders, potentially relevant after NMD inhibition.</li>
  <li><strong>Tier 3 — NMD-insensitive + IC50 &lt;500 nM:</strong> Already translated — controls for immune response without NMD inhibition.</li>
  <li><strong>Unclassified:</strong> Missense variants or insufficient transcript information.</li>
</ul>

<h2>Tier 1 — NMD-sensitive + IC50 &lt;50 nM (top 30 across cohort by IC50)</h2>
{_tbl(t1, COLS, "No Tier 1 candidates in this cohort.")}

<h2>Tier 2 — NMD-sensitive + IC50 50–500 nM (top 30 across cohort by IC50)</h2>
{_tbl(t2, COLS, "No Tier 2 candidates in this cohort.")}

<h2>Tier 3 — NMD-insensitive controls (top 30 across cohort by IC50)</h2>
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
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}\n  GBM Pipeline — Stage 3 cohort summary\n{'='*60}\n")
    print(f"[INFO] Reading per-sample dirs from: {args.input_dir}")

    n_input_samples = sum(1 for p in args.input_dir.iterdir() if p.is_dir() and SAMPLE_RE.match(p.name))
    cohort  = load_cohort(args.input_dir)
    summary = cohort_summary(cohort, n_input_samples)
    paired  = per_patient_paired(cohort)

    if not cohort.empty:
        cohort.to_csv(args.out_dir / "cohort_candidates.tsv", sep="\t", index=False)
        cohort[cohort["priority_tier"]=="TIER1"].to_csv(
            args.out_dir / "cohort_tier1.tsv", sep="\t", index=False)
    summary.to_csv(args.out_dir / "cohort_summary.tsv", sep="\t", index=False)
    paired.to_csv (args.out_dir / "cohort_paired.tsv",  sep="\t", index=False)

    generate_report(cohort, summary, paired, args.out_dir, n_input_samples)

    print(f"\n{'='*60}\n  Cohort summary done. Output: {args.out_dir.resolve()}\n{'='*60}\n")


if __name__ == "__main__":
    main()
