#!/usr/bin/env python3
"""
gbm_analysis.py - Stage 1 of the GBM NMD-Neoantigen Pipeline.

Builds per-sample somatic mutation summaries for 28 paired primary/recurrent
glioblastoma samples (56 SnpEff-annotated VCFs) and the paired statistics that
downstream stages and the cohort report consume.

STRUCTURE
---------
  1. VCF parsing            parse one VCF -> per-sample counts (parse_sample)
  2. Cohort summary         all VCFs -> sample_mutation_summary.tsv (build_summary)
  3. Paired overlap         T vs M variant sets per patient (compute_overlap)
  4. Paired statistics      Wilcoxon signed-rank + BH-FDR (compute_paired_stats)
  5. Plots                  seven PNGs at 300 dpi (make_*_plot)
  6. HTML report            self-contained stage1_report.html (build_report)

OUTPUTS (written to --out_dir)
------------------------------
  sample_mutation_summary.tsv     per-sample mutation summary (replaces
                                  summary_mutation_burden.tsv)
  paired_variant_overlap.tsv      T/M shared-vs-private counts per patient
  paired_stats.tsv                paired Wilcoxon tests + BH-FDR per metric
  all_truncating_variants.tsv     every truncating variant (FS/SG/splice/start-lost)
                                  (replaces all_fs_sg_variants.tsv)
  gene_recurrence.tsv             genes recurrently hit by truncating variants
  hla_typing_summary.tsv          HLA Class I alleles per sample (if --hla_csv given)
  plot_paired_tmb.png                     plot_paired_truncating_fraction.png
  plot_delta_tmb_waterfall.png            plot_delta_truncating_fraction_waterfall.png
  plot_mutation_class_stacked.png         plot_tmb_vs_truncating_fraction.png
  plot_variant_overlap_summary.png        plot_sbs11_boxplot.png (report section 7)
  stage1_report.html              self-contained HTML report (replaces report.html)

DEPRECATED / REMOVED
--------------------
  summary_mutation_burden.tsv, tmz_signature.tsv, report.html and all previous
  plot_*.png are deleted from --out_dir at the start of a run; their content now
  lives in the files above (TMZ data is in the sbs11_count / sbs11_pct columns of
  sample_mutation_summary.tsv).

USAGE
-----
  python gbm_analysis.py --input_dir data/vcf \\
                         --hla_csv data/hla_typing_classI.csv \\
                         --out_dir results/run_<TS>/1_gbm_analysis
  (--vcf_dir / --hla_dir are accepted as aliases for backward compatibility.)
"""

from __future__ import annotations

import argparse
import base64
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cyvcf2 import VCF
from scipy.stats import spearmanr, wilcoxon

# ── PALETTE ──────────────────────────────────────────────────────────────────
COLORS: Dict[str, str] = {
    "primary": "#377AB8",      # blue   - primary (T)
    "recurrent": "#D85A30",    # orange - recurrent (M)
    "shared": "#888888",       # grey   - shared variants
    "truncating": "#C0392B",   # red    - truncating
    "missense": "#F39C12",     # amber  - missense
    "synonymous": "#95A5A6",   # grey   - synonymous
}

# ── VCF FILENAME PARSING ─────────────────────────────────────────────────────
VCF_PATTERN = re.compile(r"^(\d+)([TM])-ensemble-annotated\.vcf(\.gz)?$")

# ── CONSEQUENCE PRIORITY ─────────────────────────────────────────────────────
# Worst-consequence ranking, high priority first. Every term below is "coding";
# anything not in this list is treated as non-coding.
CONSEQUENCE_PRIORITY: List[str] = [
    "frameshift_variant",
    "stop_gained",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
    "stop_lost",
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "splice_region_variant",
    "synonymous_variant",
]
_RANK: Dict[str, int] = {term: i for i, term in enumerate(CONSEQUENCE_PRIORITY)}
TRUNCATING_SPLICE = {"splice_donor_variant", "splice_acceptor_variant"}

# Metrics compared between paired timepoints (paired_stats.tsv).
PAIRED_METRICS: List[str] = [
    "snv", "indel", "coding", "missense",
    "truncating_total", "truncating_fraction", "tmb_per_mb", "sbs11_pct",
]

# Output files deleted at the start of a run (superseded by the new outputs).
DEPRECATED_OUTPUTS: List[str] = [
    "summary_mutation_burden.tsv", "tmz_signature.tsv", "report.html",
    "all_fs_sg_variants.tsv",
    "plot_mutation_burden.png", "plot_fs_sg_burden.png", "plot_paired_scatter.png",
    "plot_snv_indel.png", "plot_variant_overlap.png", "plot_m_only_fs_sg.png",
    "plot_tmz_signature.png", "plot_top_recurrent_genes.png",
]


# ═════════════════════════════════════════════════════════════════════════════
# FASTA RANDOM ACCESS
# ═════════════════════════════════════════════════════════════════════════════

class FastaReader:
    """Context-managed random-access reader for an indexed (.fai) FASTA.

    Inputs:  path to a bgzip-free FASTA with a sibling ``<path>.fai`` index.
    Provides ``base(chrom, pos)`` returning the single 1-based reference base
    (upper-case) at a locus, or ``None`` if out of range / contig unknown.
    Transparently reconciles UCSC-style ('chr1') VCF contigs with Ensembl-style
    ('1') FASTA contigs by trying both forms. Reused across SBS11 and any future
    signature work. Side effect: holds one open binary file handle until closed.
    """

    def __init__(self, fasta_path: Path) -> None:
        self.path = Path(fasta_path)
        fai_path = Path(str(self.path) + ".fai")
        if not self.path.exists() or not fai_path.exists():
            raise FileNotFoundError(f"FASTA or .fai index missing for {self.path}")
        self._index: Dict[str, Dict[str, int]] = {}
        with open(fai_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 5:
                    self._index[parts[0]] = {
                        "length": int(parts[1]), "offset": int(parts[2]),
                        "bases": int(parts[3]), "bytes": int(parts[4]),
                    }
        self._name_cache: Dict[str, Optional[str]] = {}
        self._fh = None

    # -- context manager ------------------------------------------------------
    def __enter__(self) -> "FastaReader":
        self._fh = open(self.path, "rb")
        return self

    def __exit__(self, *exc) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None

    # -- contig name reconciliation ------------------------------------------
    def _resolve(self, chrom: str) -> Optional[str]:
        if chrom in self._name_cache:
            return self._name_cache[chrom]
        candidates = [chrom,
                      chrom[3:] if chrom.startswith("chr") else "chr" + chrom]
        resolved = next((c for c in candidates if c in self._index), None)
        self._name_cache[chrom] = resolved
        return resolved

    # -- single-base lookup ---------------------------------------------------
    def base(self, chrom: str, pos: int) -> Optional[str]:
        """Return the upper-case reference base at 1-based ``pos`` or None."""
        name = self._resolve(chrom)
        if name is None or pos < 1:
            return None
        info = self._index[name]
        if pos > info["length"]:
            return None
        zero = pos - 1
        line, rem = divmod(zero, info["bases"])
        offset = info["offset"] + line * info["bytes"] + rem
        self._fh.seek(offset)
        ch = self._fh.read(1).decode("ascii", errors="ignore")
        return ch.upper() if ch.isalpha() else None


# ═════════════════════════════════════════════════════════════════════════════
# SnpEff ANNOTATION PARSING
# ═════════════════════════════════════════════════════════════════════════════

def worst_consequence(ann: str) -> Tuple[Optional[str], Optional[str]]:
    """Extract the worst (highest-priority) coding consequence from an ANN field.

    Input:  raw SnpEff ``ANN=`` value (comma-separated entries, pipe-delimited
            fields; field 2 is the '&'-joined consequence term list).
    Output: ``(term, gene)`` where ``term`` is the highest-priority coding
            consequence in CONSEQUENCE_PRIORITY (or None if the variant is
            non-coding), and ``gene`` is the gene symbol of the entry carrying
            that term (or None).
    """
    if not ann:
        return None, None
    best_rank: Optional[int] = None
    best_gene: Optional[str] = None
    for entry in ann.split(","):
        fields = entry.split("|")
        if len(fields) < 2:
            continue
        gene = fields[3] if len(fields) > 3 else None
        for term in fields[1].split("&"):
            rank = _RANK.get(term)
            if rank is not None and (best_rank is None or rank < best_rank):
                best_rank, best_gene = rank, gene
    if best_rank is None:
        return None, None
    return CONSEQUENCE_PRIORITY[best_rank], best_gene


def _is_pass(variant) -> bool:
    """True if FILTER is PASS / '.' / empty (cyvcf2 maps all of these to None)."""
    return (not variant.FILTER) or variant.FILTER == "PASS"


# ═════════════════════════════════════════════════════════════════════════════
# PER-SAMPLE PARSING
# ═════════════════════════════════════════════════════════════════════════════

def parse_sample(vcf_path: Path, sample: str, patient: str, timepoint: str,
                 fasta: Optional[FastaReader], tmb_denom: float
                 ) -> Tuple[Dict[str, object], pd.DataFrame, Set[str]]:
    """Parse one SnpEff-annotated VCF into a per-sample summary row.

    Inputs:  ``vcf_path`` VCF file; ``sample`` ("11_T"), ``patient`` ("11"),
             ``timepoint`` ("primary"/"recurrent") labels; ``fasta`` an open
             FastaReader for SBS11 context (or None to skip SBS11); ``tmb_denom``
             callable megabases for TMB normalisation.
    Outputs: ``(row, truncating_df, pass_keys)`` - the summary dict with all
             columns of sample_mutation_summary.tsv, a DataFrame of truncating
             variants for all_truncating_variants.tsv, and the set of PASS
             variant identity keys ("chrom:pos:ref>alt") for overlap analysis.
    Side effects: reads ``vcf_path`` (streamed once) and ``fasta`` (random access).
    """
    total = pass_n = snv = indel = 0
    coding = synonymous = missense = 0
    tr_fs = tr_sg = tr_splice = tr_start = 0
    sbs11 = 0
    keys: Set[str] = set()
    truncating_rows: List[Dict[str, object]] = []

    for v in VCF(str(vcf_path)):
        total += 1
        if not _is_pass(v):
            continue
        pass_n += 1
        ref, alt = v.REF.upper(), v.ALT[0].upper()
        keys.add(f"{v.CHROM}:{v.POS}:{ref}>{alt}")

        is_snv = len(ref) == 1 and len(alt) == 1
        if is_snv:
            snv += 1
        else:
            indel += 1

        term, gene = worst_consequence(str(v.INFO.get("ANN", "") or ""))
        if term is not None:           # term in CONSEQUENCE_PRIORITY ⇒ coding
            coding += 1
            if term == "synonymous_variant":
                synonymous += 1
            elif term == "missense_variant":
                missense += 1
            elif term == "frameshift_variant":
                tr_fs += 1
            elif term == "stop_gained":
                tr_sg += 1
            elif term in TRUNCATING_SPLICE:
                tr_splice += 1
            elif term == "start_lost":
                tr_start += 1
            if term in ("frameshift_variant", "stop_gained", "start_lost") \
                    or term in TRUNCATING_SPLICE:
                truncating_rows.append({
                    "sample": sample, "patient": patient, "timepoint": timepoint,
                    "chrom": v.CHROM, "pos": v.POS, "ref": ref, "alt": alt,
                    "gene": gene or "Unknown", "consequence": term})

        # SBS11 (TMZ): C>T with 3' ref base C, or strand-complement G>A with 5'
        # ref base G (NpCpC trinucleotide context).
        if is_snv and fasta is not None:
            if ref == "C" and alt == "T":
                if fasta.base(v.CHROM, v.POS + 1) == "C":
                    sbs11 += 1
            elif ref == "G" and alt == "A":
                if fasta.base(v.CHROM, v.POS - 1) == "G":
                    sbs11 += 1

    tr_total = tr_fs + tr_sg + tr_splice + tr_start
    row: Dict[str, object] = {
        "sample": sample, "patient": patient, "timepoint": timepoint,
        "total_variants": total, "pass_variants": pass_n,
        "snv": snv, "indel": indel,
        "coding": coding, "synonymous": synonymous, "missense": missense,
        "nonsynonymous": coding - synonymous,
        "truncating_fs": tr_fs, "truncating_sg": tr_sg,
        "truncating_splice": tr_splice, "truncating_startlost": tr_start,
        "truncating_total": tr_total,
        "truncating_fraction": round(tr_total / coding, 4) if coding else 0.0,
        "tmb_per_mb": round(pass_n / tmb_denom, 4),
        "sbs11_count": sbs11,
        "sbs11_pct": round(sbs11 / snv * 100, 1) if snv else 0.0,
    }
    truncating_df = pd.DataFrame(truncating_rows)
    return row, truncating_df, keys


# ═════════════════════════════════════════════════════════════════════════════
# COHORT SUMMARY
# ═════════════════════════════════════════════════════════════════════════════

def _list_vcfs(input_dir: Path) -> List[Path]:
    """Return sorted VCF paths (plain then gzipped) in ``input_dir``."""
    return sorted(input_dir.glob("*.vcf")) + sorted(input_dir.glob("*.vcf.gz"))


def build_summary(input_dir: Path, fasta: Optional[FastaReader], tmb_denom: float
                  ) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Set[str]]]:
    """Parse every VCF in ``input_dir`` into the cohort summary.

    Inputs:  ``input_dir`` directory of *-ensemble-annotated.vcf files; an open
             ``fasta`` (or None); ``tmb_denom`` callable megabases.
    Outputs: ``(summary_df, truncating_df, keys_by_sample)`` - the per-sample
             summary table (columns in spec order), the concatenated truncating
             variant table, and a {sample: set(variant keys)} map for overlap.
    Side effects: prints per-sample progress; exits if no VCFs are found.
    """
    print("[STEP 1] Parsing VCFs ...")
    vcfs = _list_vcfs(input_dir)
    if not vcfs:
        sys.exit(f"[ERROR] No VCF files found in {input_dir}")
    print(f"  {len(vcfs)} VCF file(s) found")

    rows: List[Dict[str, object]] = []
    truncating_frames: List[pd.DataFrame] = []
    keys_by_sample: Dict[str, Set[str]] = {}

    for vcf_path in vcfs:
        m = VCF_PATTERN.match(vcf_path.name)
        if not m:
            print(f"  [skip] {vcf_path.name} (unrecognised name)")
            continue
        patient, tp = m.group(1), m.group(2)
        sample = f"{patient}_{tp}"
        timepoint = "primary" if tp == "T" else "recurrent"
        row, tr_df, keys = parse_sample(vcf_path, sample, patient, timepoint,
                                        fasta, tmb_denom)
        rows.append(row)
        keys_by_sample[sample] = keys
        if not tr_df.empty:
            truncating_frames.append(tr_df)
        print(f"  {sample:>6}  PASS={row['pass_variants']:>6}  "
              f"trunc={row['truncating_total']:>4}  SBS11%={row['sbs11_pct']}")

    summary_cols = [
        "sample", "patient", "timepoint", "total_variants", "pass_variants",
        "snv", "indel", "coding", "synonymous", "missense", "nonsynonymous",
        "truncating_fs", "truncating_sg", "truncating_splice",
        "truncating_startlost", "truncating_total", "truncating_fraction",
        "tmb_per_mb", "sbs11_count", "sbs11_pct",
    ]
    summary = pd.DataFrame(rows, columns=summary_cols)
    summary = summary.sort_values(
        ["patient", "timepoint"],
        key=lambda s: s.map(int) if s.name == "patient"
        else s.map({"primary": 0, "recurrent": 1})).reset_index(drop=True)

    truncating = (pd.concat(truncating_frames, ignore_index=True)
                  if truncating_frames else pd.DataFrame(
                      columns=["sample", "patient", "timepoint", "chrom", "pos",
                               "ref", "alt", "gene", "consequence"]))
    return summary, truncating, keys_by_sample


# ═════════════════════════════════════════════════════════════════════════════
# PAIRED VARIANT OVERLAP
# ═════════════════════════════════════════════════════════════════════════════

def compute_overlap(keys_by_sample: Dict[str, Set[str]]) -> pd.DataFrame:
    """Compute T-vs-M variant-set overlap for patients with both timepoints.

    Input:  {sample ("11_T"/"11_M"): set of "chrom:pos:ref>alt" PASS keys}.
    Output: DataFrame with columns patient, t_total, m_total, shared, t_private,
            m_private, shared_pct_of_t, shared_pct_of_m. Singleton patients
            (only one timepoint) are skipped.
    """
    print("[STEP 3] Computing paired variant overlap ...")
    patients: Dict[str, Dict[str, Set[str]]] = {}
    for sample, keys in keys_by_sample.items():
        patient, tp = sample.split("_")
        patients.setdefault(patient, {})[tp] = keys

    rows: List[Dict[str, object]] = []
    for patient in sorted(patients, key=int):
        tp_map = patients[patient]
        if "T" not in tp_map or "M" not in tp_map:
            continue
        t_keys, m_keys = tp_map["T"], tp_map["M"]
        shared = len(t_keys & m_keys)
        rows.append({
            "patient": patient,
            "t_total": len(t_keys), "m_total": len(m_keys), "shared": shared,
            "t_private": len(t_keys - m_keys), "m_private": len(m_keys - t_keys),
            "shared_pct_of_t": round(shared / len(t_keys) * 100, 1) if t_keys else 0.0,
            "shared_pct_of_m": round(shared / len(m_keys) * 100, 1) if m_keys else 0.0,
        })
    df = pd.DataFrame(rows)
    print(f"  {len(df)} paired patient(s)")
    return df


# ═════════════════════════════════════════════════════════════════════════════
# PAIRED STATISTICS
# ═════════════════════════════════════════════════════════════════════════════

def _bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment; NaN p-values are excluded and kept NaN."""
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


def compute_paired_stats(summary: pd.DataFrame) -> pd.DataFrame:
    """Paired Wilcoxon signed-rank tests across timepoints for each metric.

    Input:  the per-sample summary table.
    Output: DataFrame (one row per metric in PAIRED_METRICS) with columns metric,
            n_pairs, median_primary, median_recurrent, median_delta, mean_delta,
            wilcoxon_stat, p_value, p_value_adj_bh, note. The BH-FDR correction is
            applied across all rows of this single table (one correction family).
            p_value is NaN (with a note) when a metric has <6 pairs or zero
            variance in the paired deltas.
    """
    print("[STEP 4] Computing paired statistics (Wilcoxon + BH-FDR) ...")
    prim = summary[summary.timepoint == "primary"].set_index("patient")
    recu = summary[summary.timepoint == "recurrent"].set_index("patient")
    paired_patients = sorted(prim.index.intersection(recu.index), key=int)

    rows: List[Dict[str, object]] = []
    for metric in PAIRED_METRICS:
        p_vals = prim.loc[paired_patients, metric].astype(float).to_numpy()
        r_vals = recu.loc[paired_patients, metric].astype(float).to_numpy()
        delta = r_vals - p_vals
        n_pairs = len(delta)

        stat: float = np.nan
        p_value: float = np.nan
        note = ""
        if n_pairs < 6:
            note = "fewer than 6 pairs"
        elif np.allclose(delta, delta[0]):  # zero variance (often all-zero deltas)
            note = "zero variance in deltas"
        else:
            try:
                res = wilcoxon(r_vals, p_vals)
                stat, p_value = float(res.statistic), float(res.pvalue)
            except ValueError as exc:       # e.g. all differences zero after drop
                note = f"wilcoxon undefined ({exc})"

        rows.append({
            "metric": metric, "n_pairs": n_pairs,
            "median_primary": round(float(np.median(p_vals)), 4),
            "median_recurrent": round(float(np.median(r_vals)), 4),
            "median_delta": round(float(np.median(delta)), 4),
            "mean_delta": round(float(np.mean(delta)), 4),
            "wilcoxon_stat": stat, "p_value": p_value,
            "p_value_adj_bh": np.nan, "note": note,
        })

    stats = pd.DataFrame(rows)
    stats["p_value_adj_bh"] = _bh_adjust(stats["p_value"].to_numpy(dtype=float))
    print(f"  {len(stats)} metric(s) tested over {len(paired_patients)} pairs")
    return stats


# ═════════════════════════════════════════════════════════════════════════════
# GENE-LEVEL RECURRENCE  (kept for backward compatibility; truncating definition)
# ═════════════════════════════════════════════════════════════════════════════

def analyse_gene_recurrence(truncating: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """Count patients carrying a truncating hit per gene, split by timepoint.

    Input:  the concatenated truncating-variant table (from build_summary).
    Output: DataFrame with gene, n_primary, n_recurrent, n_either, recurrent_only;
            also written to gene_recurrence.tsv.
    """
    print("[STEP 2] Gene-level recurrence ...")
    if truncating.empty:
        df = pd.DataFrame(columns=["gene", "n_primary", "n_recurrent",
                                   "n_either", "recurrent_only"])
        df.to_csv(out_dir / "gene_recurrence.tsv", sep="\t", index=False)
        return df
    hits: Dict[str, Dict[str, Set[str]]] = {}
    for _, r in truncating.iterrows():
        tp = "T" if r["timepoint"] == "primary" else "M"
        hits.setdefault(r["gene"], {"T": set(), "M": set()})[tp].add(r["patient"])
    df = pd.DataFrame([
        {"gene": g, "n_primary": len(h["T"]), "n_recurrent": len(h["M"]),
         "n_either": len(h["T"] | h["M"]), "recurrent_only": len(h["M"] - h["T"])}
        for g, h in hits.items()
    ]).sort_values(["n_recurrent", "n_either"], ascending=False).reset_index(drop=True)
    df.to_csv(out_dir / "gene_recurrence.tsv", sep="\t", index=False)
    return df


# ═════════════════════════════════════════════════════════════════════════════
# HLA TYPING  (kept unchanged in output: hla_typing_summary.tsv)
# ═════════════════════════════════════════════════════════════════════════════

def load_hla_typing(hla_path: Path) -> pd.DataFrame:
    """Load HLA Class I typing from a CSV file or a directory of Dragen TSVs.

    Input:  ``hla_path`` either a CSV (Sample_ID, HLA_Types) or a directory
            containing such a CSV or per-sample ``S<id>_<T|M>.hla.tsv`` files.
    Output: DataFrame with Sample_ID and HLA_Types columns.
    """
    if hla_path.is_file() and hla_path.suffix == ".csv":
        df = pd.read_csv(hla_path)
        print(f"  HLA CSV loaded: {len(df)} samples from {hla_path.name}")
        return df
    if hla_path.is_dir():
        csvs = list(hla_path.glob("*.csv"))
        if csvs:
            df = pd.read_csv(csvs[0])
            print(f"  HLA CSV loaded: {len(df)} samples from {csvs[0].name}")
            return df
        tsv_pat = re.compile(r"^S?(\d+)_(T|M)\.hla\.tsv$")
        records: List[Dict[str, str]] = []
        for f in sorted(hla_path.iterdir()):
            m = tsv_pat.match(f.name)
            if not m:
                continue
            sample = f"{m.group(1)}_{m.group(2)}"
            tdf = pd.read_csv(f, sep="\t")
            alleles: List[str] = []
            for gene in ("A", "B", "C"):
                grow = tdf[tdf["gene"] == gene] if "gene" in tdf else tdf.iloc[0:0]
                if grow.empty:
                    continue
                for col in ("allele_1", "allele_2"):
                    a = grow.iloc[0].get(col)
                    if pd.notna(a) and str(a) not in ("NA", ""):
                        alleles.append(f"HLA-{gene}*{a}" if "*" not in str(a)
                                       else f"HLA-{a}")
            records.append({"Sample_ID": sample, "HLA_Types": ";".join(alleles)})
        return pd.DataFrame(records)
    raise FileNotFoundError(f"HLA path not found: {hla_path}")


# ═════════════════════════════════════════════════════════════════════════════
# PLOTS
# ═════════════════════════════════════════════════════════════════════════════

def _style_axes(ax: plt.Axes) -> None:
    """Apply the clean shared plot style (drop top/right spines)."""
    ax.spines[["top", "right"]].set_visible(False)


def _save(fig: plt.Figure, path: Path) -> None:
    """Save a figure at 300 dpi and close it."""
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print(f"  [plot] {path.name}")


def _pivot(summary: pd.DataFrame, metric: str) -> Tuple[List[str], List[float], List[float]]:
    """Return (paired_patients, primary_values, recurrent_values) for ``metric``."""
    prim = summary[summary.timepoint == "primary"].set_index("patient")[metric]
    recu = summary[summary.timepoint == "recurrent"].set_index("patient")[metric]
    patients = sorted(prim.index.intersection(recu.index), key=int)
    return patients, [float(prim[p]) for p in patients], [float(recu[p]) for p in patients]


def _p_label(stats: pd.DataFrame, metric: str) -> str:
    """Format a 'p = ...' annotation for ``metric`` from paired_stats."""
    row = stats[stats.metric == metric]
    if row.empty or pd.isna(row.iloc[0]["p_value"]):
        return "p = n/a"
    p = float(row.iloc[0]["p_value"])
    return f"p = {p:.1e}" if p < 1e-3 else f"p = {p:.3f}"


def make_paired_connected_plot(summary: pd.DataFrame, stats: pd.DataFrame,
                               metric: str, ylabel: str, title: str,
                               out_path: Path) -> None:
    """Connected dot plot joining each patient's primary and recurrent value."""
    patients, prim, recu = _pivot(summary, metric)
    fig, ax = plt.subplots(figsize=(6, 7))
    for yp, yr in zip(prim, recu):
        ax.plot([0, 1], [yp, yr], color="#cccccc", lw=0.9, zorder=1)
    ax.scatter([0] * len(prim), prim, color=COLORS["primary"], s=45, zorder=3)
    ax.scatter([1] * len(recu), recu, color=COLORS["recurrent"], s=45, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Primary", "Recurrent"])
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.text(0.98, 0.98, _p_label(stats, metric), transform=ax.transAxes,
            ha="right", va="top", fontsize=11,
            bbox=dict(boxstyle="round", fc="white", ec="#cccccc"))
    _style_axes(ax)
    _save(fig, out_path)


def make_waterfall_plot(summary: pd.DataFrame, metric: str, ylabel: str,
                        title: str, out_path: Path) -> None:
    """Waterfall bar plot of (recurrent - primary) per patient, sorted by delta."""
    patients, prim, recu = _pivot(summary, metric)
    deltas = [(p, r - t) for p, t, r in zip(patients, prim, recu)]
    deltas.sort(key=lambda x: x[1])
    labels = [d[0] for d in deltas]
    values = [d[1] for d in deltas]
    colors = [COLORS["truncating"] if v > 0 else COLORS["primary"] for v in values]
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar(range(len(values)), values, color=colors, width=0.7)
    ax.axhline(0, color="#333333", lw=0.8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlabel("Patient (sorted by Δ)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    _style_axes(ax)
    _save(fig, out_path)


def make_mutation_class_stacked_plot(summary: pd.DataFrame, out_path: Path) -> None:
    """Stacked bar of synonymous / missense / truncating counts per sample."""
    df = summary  # already sorted: patients numeric, primary before recurrent
    x = range(len(df))
    syn = df["synonymous"].to_numpy()
    mis = df["missense"].to_numpy()
    tru = df["truncating_total"].to_numpy()
    fig, ax = plt.subplots(figsize=(16, 6))
    ax.bar(x, syn, color=COLORS["synonymous"], label="Synonymous")
    ax.bar(x, mis, bottom=syn, color=COLORS["missense"], label="Missense")
    ax.bar(x, tru, bottom=syn + mis, color=COLORS["truncating"], label="Truncating")
    ax.set_xticks(list(x))
    ax.set_xticklabels(df["sample"], rotation=90, fontsize=8)
    ax.set_ylabel("Coding variant count")
    ax.set_title("Mutation class composition per sample")
    ax.legend()
    _style_axes(ax)
    _save(fig, out_path)


def make_overlap_plot(overlap: pd.DataFrame, out_path: Path) -> None:
    """Grouped bars (t_private / shared / m_private) per paired patient."""
    if overlap.empty:
        return
    patients = overlap["patient"].tolist()
    x = np.arange(len(patients))
    w = 0.27
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.bar(x - w, overlap["t_private"], w, color=COLORS["primary"], label="Primary-private")
    ax.bar(x, overlap["shared"], w, color=COLORS["shared"], label="Shared")
    ax.bar(x + w, overlap["m_private"], w, color=COLORS["recurrent"], label="Recurrent-private")
    ax.set_xticks(x)
    ax.set_xticklabels(patients, fontsize=9)
    ax.set_xlabel("Patient")
    ax.set_ylabel("PASS variant count")
    ax.set_title("Variant overlap - private vs shared per patient")
    ax.legend()
    _style_axes(ax)
    _save(fig, out_path)


def make_tmb_vs_truncating_plot(summary: pd.DataFrame, out_path: Path) -> None:
    """Scatter of TMB (log x) vs truncating fraction for all samples + Spearman."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for tp, color in (("primary", COLORS["primary"]), ("recurrent", COLORS["recurrent"])):
        sub = summary[summary.timepoint == tp]
        ax.scatter(sub["tmb_per_mb"], sub["truncating_fraction"],
                   color=color, s=45, alpha=0.8, label=tp.capitalize())
    x = summary["tmb_per_mb"].to_numpy(dtype=float)
    y = summary["truncating_fraction"].to_numpy(dtype=float)
    valid = x > 0
    if valid.sum() > 2:
        rho, p = spearmanr(x[valid], y[valid])
        ann = f"Spearman ρ = {rho:.2f}\np = {p:.1e}" if p < 1e-3 else f"Spearman ρ = {rho:.2f}\np = {p:.3f}"
        ax.text(0.98, 0.98, ann, transform=ax.transAxes, ha="right", va="top",
                fontsize=11, bbox=dict(boxstyle="round", fc="white", ec="#cccccc"))
    ax.set_xscale("log")
    ax.set_xlabel("TMB per Mb (log scale)")
    ax.set_ylabel("Truncating fraction")
    ax.set_title("TMB vs truncating fraction (all samples)")
    ax.legend()
    _style_axes(ax)
    _save(fig, out_path)


def make_sbs11_boxplot(summary: pd.DataFrame, out_path: Path) -> None:
    """Boxplot of sbs11_pct split by timepoint (report section 7)."""
    prim = summary[summary.timepoint == "primary"]["sbs11_pct"].to_numpy(dtype=float)
    recu = summary[summary.timepoint == "recurrent"]["sbs11_pct"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(6, 6))
    bp = ax.boxplot([prim, recu], tick_labels=["Primary", "Recurrent"],
                    patch_artist=True, widths=0.5)
    for patch, color in zip(bp["boxes"], (COLORS["primary"], COLORS["recurrent"])):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    for med in bp["medians"]:
        med.set_color("#222222")
    # jittered points
    for i, vals in enumerate((prim, recu), start=1):
        jitter = np.linspace(-0.08, 0.08, len(vals)) if len(vals) else []
        ax.scatter(np.full(len(vals), i) + jitter, vals, color="#333333",
                   s=18, alpha=0.6, zorder=3)
    ax.set_ylabel("SBS11 % of SNVs")
    ax.set_title("TMZ signature (SBS11) by timepoint")
    _style_axes(ax)
    _save(fig, out_path)


def make_all_plots(summary: pd.DataFrame, overlap: pd.DataFrame,
                   stats: pd.DataFrame, out_dir: Path) -> None:
    """Generate the full Stage 1 plot set (300 dpi PNGs) into ``out_dir``."""
    print("[STEP 5] Generating plots ...")
    make_paired_connected_plot(summary, stats, "tmb_per_mb", "TMB per Mb",
                               "Paired TMB: primary vs recurrent",
                               out_dir / "plot_paired_tmb.png")
    make_waterfall_plot(summary, "tmb_per_mb", "Δ TMB per Mb (recurrent − primary)",
                        "TMB change at recurrence",
                        out_dir / "plot_delta_tmb_waterfall.png")
    make_mutation_class_stacked_plot(summary, out_dir / "plot_mutation_class_stacked.png")
    make_overlap_plot(overlap, out_dir / "plot_variant_overlap_summary.png")
    make_paired_connected_plot(summary, stats, "truncating_fraction",
                               "Truncating fraction",
                               "Paired truncating fraction: primary vs recurrent",
                               out_dir / "plot_paired_truncating_fraction.png")
    make_waterfall_plot(summary, "truncating_fraction",
                        "Δ truncating fraction (recurrent − primary)",
                        "Truncating fraction change at recurrence",
                        out_dir / "plot_delta_truncating_fraction_waterfall.png")
    make_tmb_vs_truncating_plot(summary, out_dir / "plot_tmb_vs_truncating_fraction.png")
    make_sbs11_boxplot(summary, out_dir / "plot_sbs11_boxplot.png")


# ═════════════════════════════════════════════════════════════════════════════
# HTML REPORT
# ═════════════════════════════════════════════════════════════════════════════

_CSS = """
* { box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
       color: #222; background: #fafafa; margin: 0; font-size: 13px; line-height: 1.5; }
.wrap { max-width: 1040px; margin: 0 auto; padding: 32px 24px 64px; }
header { border-bottom: 2px solid #e3e3e3; padding-bottom: 18px; margin-bottom: 28px; }
header h1 { font-size: 23px; margin: 0 0 6px; font-weight: 600; }
header p { color: #777; margin: 2px 0; font-size: 12px; }
.note { background: #fff8e6; border: 1px solid #f0e2b6; border-radius: 6px;
        padding: 10px 14px; font-size: 12px; color: #6b5a1e; margin-top: 12px; }
h2 { font-size: 16px; font-weight: 600; border-bottom: 1px solid #e5e5e5;
     padding-bottom: 6px; margin: 36px 0 14px; }
.cards { display: flex; gap: 12px; flex-wrap: wrap; margin: 14px 0 4px; }
.card { background: #fff; border: 1px solid #e5e5e5; border-radius: 8px;
        padding: 12px 18px; min-width: 130px; }
.card .lbl { font-size: 11px; color: #888; }
.card .val { font-size: 22px; font-weight: 600; }
.tablewrap { max-height: 460px; overflow: auto; border: 1px solid #e5e5e5;
             border-radius: 8px; margin: 12px 0; }
table { width: 100%; border-collapse: collapse; font-size: 12px; }
th { position: sticky; top: 0; background: #f1f3f5; text-align: right;
     padding: 7px 10px; font-weight: 600; border-bottom: 2px solid #dcdcdc; z-index: 1; }
th:first-child, td:first-child { text-align: left; }
td { padding: 5px 10px; border-bottom: 1px solid #eee; text-align: right; }
tbody tr:nth-child(even) { background: #fafbfc; }
img { max-width: 100%; height: auto; display: block; margin: 14px auto;
      border-radius: 6px; box-shadow: 0 1px 4px rgba(0,0,0,.10); }
.cap { font-size: 11px; color: #999; text-align: center; margin: -6px 0 18px; }
"""


def _b64_img(path: Path) -> str:
    """Return a base64 data URI for ``path`` or '' if it does not exist."""
    if not path.exists():
        return ""
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()


def _fig(out_dir: Path, name: str, caption: str) -> str:
    """Build an embedded <img> + caption block for plot ``name``."""
    uri = _b64_img(out_dir / name)
    if not uri:
        return f"<p class='cap'>[missing plot: {name}]</p>"
    return f'<img src="{uri}" alt="{name}"><p class="cap">{caption}</p>'


def _table(df: pd.DataFrame, float_cols_1dp: Optional[Set[str]] = None) -> str:
    """Render a DataFrame as a scrollable HTML table with a sticky header."""
    float_cols_1dp = float_cols_1dp or set()

    def fmt(col: str, v: object) -> str:
        if isinstance(v, float):
            if np.isnan(v):
                return "-"
            if col in float_cols_1dp:
                return f"{v:.1f}"
            return f"{v:.4g}" if v != int(v) else f"{int(v)}"
        return str(v)

    head = "".join(f"<th>{c}</th>" for c in df.columns)
    body = "".join(
        "<tr>" + "".join(f"<td>{fmt(c, r[c])}</td>" for c in df.columns) + "</tr>"
        for _, r in df.iterrows())
    return f'<div class="tablewrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def build_report(out_dir: Path, summary: pd.DataFrame, overlap: pd.DataFrame,
                 stats: pd.DataFrame, hla_df: Optional[pd.DataFrame],
                 run_ts: str, tmb_denom: float) -> Path:
    """Assemble the self-contained Stage 1 HTML report (stage1_report.html).

    Inputs:  all Stage 1 tables plus run metadata. Images are read back from
             ``out_dir`` and base64-embedded; CSS is inlined.
    Output:  path to the written stage1_report.html.
    Side effect: writes the HTML file.
    """
    print("[STEP 6] Building HTML report ...")
    n_samples = len(summary)
    n_patients = summary["patient"].nunique()
    n_pairs = int(overlap.shape[0])
    n_primary = int((summary.timepoint == "primary").sum())
    n_recurrent = int((summary.timepoint == "recurrent").sum())

    def cards(*items: Tuple[str, object]) -> str:
        return '<div class="cards">' + "".join(
            f'<div class="card"><div class="lbl">{l}</div><div class="val">{v}</div></div>'
            for l, v in items) + "</div>"

    pct1 = {"truncating_fraction", "sbs11_pct", "shared_pct_of_t", "shared_pct_of_m"}

    sections: List[str] = []

    # Section 1 - cohort overview
    sections.append(
        "<h2>1. Cohort overview</h2>"
        + cards(("Samples", n_samples), ("Patients", n_patients),
                ("Paired patients", n_pairs), ("Primary (T)", n_primary),
                ("Recurrent (M)", n_recurrent)))

    # Section 2 - per-sample summary table
    sections.append(
        "<h2>2. Per-sample mutation summary</h2>"
        "<p>Full per-sample table (scrollable). TMB normalised to "
        f"{tmb_denom:g} Mb.</p>"
        + _table(summary, pct1))

    # Section 3 - paired statistics
    sections.append(
        "<h2>3. Paired statistics</h2>"
        "<p>Paired Wilcoxon signed-rank tests across timepoints; "
        "p-values BH-FDR corrected across the rows of this table.</p>"
        + _table(stats, {"truncating_fraction"}))

    # Section 4 - mutation burden
    sections.append(
        "<h2>4. Mutation burden</h2>"
        + _fig(out_dir, "plot_paired_tmb.png",
               "Fig 1. Paired TMB per Mb (one line per patient).")
        + _fig(out_dir, "plot_delta_tmb_waterfall.png",
               "Fig 2. ΔTMB (recurrent − primary), sorted.")
        + _fig(out_dir, "plot_mutation_class_stacked.png",
               "Fig 3. Mutation-class composition per sample."))

    # Section 5 - variant overlap
    sections.append(
        "<h2>5. Variant overlap</h2>"
        + _fig(out_dir, "plot_variant_overlap_summary.png",
               "Fig 4. Private vs shared variants per paired patient.")
        + _table(overlap, pct1))

    # Section 6 - truncating enrichment
    sections.append(
        "<h2>6. Truncating mutation enrichment</h2>"
        + _fig(out_dir, "plot_paired_truncating_fraction.png",
               "Fig 5. Paired truncating fraction.")
        + _fig(out_dir, "plot_delta_truncating_fraction_waterfall.png",
               "Fig 6. Δ truncating fraction (recurrent − primary), sorted.")
        + _fig(out_dir, "plot_tmb_vs_truncating_fraction.png",
               "Fig 7. TMB (log) vs truncating fraction, all 56 samples."))

    # Section 7 - TMZ signature
    sbs = summary[["sample", "timepoint", "sbs11_count", "sbs11_pct"]].copy()
    med_p = summary[summary.timepoint == "primary"]["sbs11_pct"].median()
    med_r = summary[summary.timepoint == "recurrent"]["sbs11_pct"].median()
    sections.append(
        "<h2>7. TMZ signature (SBS11)</h2>"
        "<p>SBS11 (temozolomide) C&gt;T in NpCpC trinucleotide context, as a "
        "percentage of SNVs. Median primary "
        f"<strong>{med_p:.1f}%</strong> vs recurrent <strong>{med_r:.1f}%</strong>.</p>"
        + _fig(out_dir, "plot_sbs11_boxplot.png",
               "Fig 8. SBS11 % of SNVs by timepoint.")
        + _table(sbs, {"sbs11_pct"}))

    hla_line = (f"{len(hla_df)} sample(s) HLA-typed" if hla_df is not None
                else "HLA typing not provided")
    header = (
        '<header><h1>GBM NMD-Neoantigen Pipeline - Stage 1</h1>'
        f'<p>Somatic mutation landscape · run {run_ts}</p>'
        f'<p>{n_samples} samples · {n_patients} patients · {n_pairs} paired · {hla_line}</p>'
        f'<div class="note">TMB normalised to a default <strong>{tmb_denom:g} Mb</strong> '
        'callable region (standard whole-exome assumption). When the capture-kit '
        'BED becomes available, set <code>TMB_DENOMINATOR_MB</code> in config.sh to '
        'the exact callable size and re-run; all TMB values scale accordingly.</div>'
        '</header>')

    html = (f'<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>GBM Stage 1 - {run_ts}</title><style>{_CSS}</style></head>'
            f'<body><div class="wrap">{header}{"".join(sections)}</div></body></html>')

    out = out_dir / "stage1_report.html"
    out.write_text(html, encoding="utf-8")
    print(f"  [report] {out.name}")
    return out


# ═════════════════════════════════════════════════════════════════════════════
# OUTPUT HOUSEKEEPING
# ═════════════════════════════════════════════════════════════════════════════

def remove_deprecated(out_dir: Path) -> None:
    """Delete superseded Stage 1 outputs (old summary, TMZ table, report, plots)."""
    removed = 0
    for name in DEPRECATED_OUTPUTS:
        p = out_dir / name
        if p.exists():
            p.unlink()
            removed += 1
    if removed:
        print(f"[INFO] Removed {removed} deprecated output file(s)")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main() -> None:
    """CLI entry point: parse VCFs, build tables/plots/report into --out_dir."""
    parser = argparse.ArgumentParser(
        description="GBM NMD-Neoantigen Pipeline - Stage 1 (mutation landscape)")
    parser.add_argument("--input_dir", "--vcf_dir", dest="input_dir", required=True,
                        help="Directory of *-ensemble-annotated.vcf files")
    parser.add_argument("--out_dir", default="./gbm_output", help="Output directory")
    parser.add_argument("--hla_csv", "--hla_dir", dest="hla", default=None,
                        help="HLA typing CSV file or directory (optional)")
    parser.add_argument("--fasta", default=None,
                        help="Indexed GRCh38 FASTA for SBS11 (default: "
                             "<input_dir>/../reference/GRCh38.primary_assembly.fa)")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = (Path(args.fasta) if args.fasta
                  else input_dir.parent / "reference" / "GRCh38.primary_assembly.fa")
    tmb_denom = float(os.environ.get("TMB_DENOMINATOR_MB", "30"))

    # Derive run timestamp from a run_<TS> output dir if possible.
    run_ts = next((part.replace("run_", "") for part in reversed(out_dir.parts)
                   if part.startswith("run_")), datetime.now().strftime("%Y%m%d_%H%M%S"))

    print("=" * 64)
    print("  GBM NMD-Neoantigen Pipeline - Stage 1")
    print(f"  input  : {input_dir}")
    print(f"  output : {out_dir}")
    print(f"  TMB denominator: {tmb_denom:g} Mb")
    print("=" * 64)

    remove_deprecated(out_dir)

    # Open the FASTA once and reuse for all SBS11 context lookups.
    fasta_cm = None
    try:
        fasta_cm = FastaReader(fasta_path).__enter__()
        print(f"[INFO] FASTA opened for SBS11 context: {fasta_path.name}")
    except FileNotFoundError:
        print(f"[WARN] FASTA/index not found at {fasta_path}; SBS11 counts = 0")

    try:
        summary, truncating, keys_by_sample = build_summary(input_dir, fasta_cm, tmb_denom)
    finally:
        if fasta_cm is not None:
            fasta_cm.__exit__(None, None, None)

    summary.to_csv(out_dir / "sample_mutation_summary.tsv", sep="\t", index=False)
    truncating.to_csv(out_dir / "all_truncating_variants.tsv", sep="\t", index=False)

    analyse_gene_recurrence(truncating, out_dir)

    overlap = compute_overlap(keys_by_sample)
    overlap.to_csv(out_dir / "paired_variant_overlap.tsv", sep="\t", index=False)

    stats = compute_paired_stats(summary)
    stats.to_csv(out_dir / "paired_stats.tsv", sep="\t", index=False)

    make_all_plots(summary, overlap, stats, out_dir)

    hla_df: Optional[pd.DataFrame] = None
    if args.hla:
        print("[STEP 7] Loading HLA typing ...")
        try:
            hla_df = load_hla_typing(Path(args.hla))
            hla_df.to_csv(out_dir / "hla_typing_summary.tsv", sep="\t", index=False)
        except Exception as exc:                       # noqa: BLE001 - non-fatal
            print(f"  [WARN] HLA loading failed: {exc}")

    build_report(out_dir, summary, overlap, stats, hla_df, run_ts, tmb_denom)

    n_pairs = int(overlap.shape[0])
    print(f"\n[OK] Stage 1 complete: {len(summary)} samples processed, "
          f"{n_pairs} paired patients, see {out_dir}/")


if __name__ == "__main__":
    main()
