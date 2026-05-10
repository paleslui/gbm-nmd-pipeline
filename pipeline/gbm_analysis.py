#!/usr/bin/env python3
"""
gbm_analysis.py — GBM Somatic Mutation Landscape Analysis (Stages 1-4)
Part of the GBM NMD-Neoantigen Pipeline.
https://github.com/paleslui/gbm-nmd-pipeline

Usage:
  python gbm_analysis.py --vcf_dir <path> --out_dir <path> [--hla_dir <path>] [--fasta <path>]
"""

import argparse, re, sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from cyvcf2 import VCF

# ── CONSTANTS ─────────────────────────────────────────────────────────────────
COL_T   = "#378ADD"   # blue  — primary
COL_M   = "#D85A30"   # coral — recurrent
COL_HI  = "#F0C040"   # amber — outlier
HIGH_IMPACT = re.compile(r"frameshift|stop_gained|stop_lost|start_lost|splice_acceptor|splice_donor")
VCF_PATTERN = re.compile(r"^(\d+)(T|M)-ensemble-annotated\.vcf(\.gz)?$")
HYPERMUTATOR_THRESHOLD = 5000


# ═════════════════════════════════════════════════════════════════════════════
# 1. VCF PARSING
#    build_summary: parse all VCFs → per-sample count DataFrame
#    compare_paired_variants: T vs M per patient → overlap DataFrame
# ═════════════════════════════════════════════════════════════════════════════

def _list_vcfs(vcf_dir: Path) -> list:
    return sorted(vcf_dir.glob("*.vcf")) + sorted(vcf_dir.glob("*.vcf.gz"))


def parse_vcf(vcf_path: Path) -> dict:
    """Parse a single SnpEff-annotated VCF. Returns count dict."""
    total = pass_total = pass_snv = pass_indel = pass_fs_sg = 0
    variants = []
    for v in VCF(str(vcf_path)):
        total += 1
        if v.FILTER and v.FILTER != "PASS":
            continue
        ann = str(v.INFO.get("ANN", ""))
        pass_total += 1
        is_indel = len(v.REF) != len(v.ALT[0])
        pass_indel += is_indel
        pass_snv   += not is_indel
        if HIGH_IMPACT.search(ann):
            pass_fs_sg += 1
            variants.append({"chrom": v.CHROM, "pos": v.POS,
                              "ref": v.REF, "alt": v.ALT[0], "ann": ann[:200]})
    return {"total": total, "pass_total": pass_total, "pass_snv": pass_snv,
            "pass_indel": pass_indel, "pass_fs_sg": pass_fs_sg, "variants": variants}


def build_summary(vcf_dir: Path) -> tuple:
    """Parse all VCFs → (summary_df, fs_sg_dict)."""
    vcfs = _list_vcfs(vcf_dir)
    if not vcfs:
        sys.exit(f"[ERROR] No VCF files in {vcf_dir}")
    print(f"[INFO] {len(vcfs)} VCF files found")
    records, fs_sg_dict = [], {}
    for vcf_path in vcfs:
        m = VCF_PATTERN.match(vcf_path.name)
        if not m:
            continue
        pid, tp, label = m.group(1), m.group(2), f"{m.group(1)}{m.group(2)}"
        timepoint = "primary" if tp == "T" else "recurrent"
        print(f"  {label} ...", end=" ", flush=True)
        s = parse_vcf(vcf_path)
        print(f"PASS={s['pass_total']}  FS/SG={s['pass_fs_sg']}")
        records.append({"sample": label, "patient": pid, "timepoint": timepoint,
                        "total": s["total"], "pass_total": s["pass_total"],
                        "pass_snv": s["pass_snv"], "pass_indel": s["pass_indel"],
                        "pass_fs_sg": s["pass_fs_sg"]})
        if s["variants"]:
            df = pd.DataFrame(s["variants"])
            df.insert(0, "sample", label); df.insert(1, "patient", pid); df.insert(2, "timepoint", timepoint)
            fs_sg_dict[label] = df
    return pd.DataFrame(records), fs_sg_dict


def compare_paired_variants(vcf_dir: Path) -> pd.DataFrame:
    """For each patient, compare T vs M variant sets → overlap DataFrame."""
    pairs = {}
    for vcf_path in _list_vcfs(vcf_dir):
        m = VCF_PATTERN.match(vcf_path.name)
        if m:
            pairs.setdefault(m.group(1), {})[m.group(2)] = vcf_path
    records = []
    for pid in sorted(pairs, key=int):
        p = pairs[pid]
        if "T" not in p or "M" not in p:
            continue
        def keys(path):
            return {f"{v.CHROM}:{v.POS}:{v.REF}>{v.ALT[0]}"
                    for v in VCF(str(path)) if not v.FILTER or v.FILTER == "PASS"}
        t_keys, m_keys = keys(p["T"]), keys(p["M"])
        m_only = m_keys - t_keys
        fs_sg = sum(1 for v in VCF(str(p["M"]))
                    if (not v.FILTER or v.FILTER == "PASS")
                    and f"{v.CHROM}:{v.POS}:{v.REF}>{v.ALT[0]}" in m_only
                    and HIGH_IMPACT.search(str(v.INFO.get("ANN", ""))))
        pct = round(fs_sg / len(m_only) * 100, 1) if m_only else 0.0
        records.append({"patient": pid, "t_total": len(t_keys), "m_total": len(m_keys),
                        "shared": len(t_keys & m_keys), "t_only": len(t_keys - m_keys),
                        "m_only": len(m_only), "m_only_fs_sg": fs_sg, "m_only_fs_sg_pct": pct})
        print(f"  Patient {pid}: M-only={len(m_only)}  FS/SG={fs_sg} ({pct}%)")
    return pd.DataFrame(records)


# ═════════════════════════════════════════════════════════════════════════════
# 2. MUTATION BURDEN — total PASS variants per sample
# ═════════════════════════════════════════════════════════════════════════════

def plot_mutation_burden(summary: pd.DataFrame, out_dir: Path):
    prim = summary[summary.timepoint == "primary"].set_index("patient")["pass_total"]
    recu = summary[summary.timepoint == "recurrent"].set_index("patient")["pass_total"]
    patients = sorted(prim.index.union(recu.index), key=int)
    x, w = range(len(patients)), 0.4
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar([i-w/2 for i in x], [prim.get(p, 0) for p in patients], w, color=COL_T, label="Primary (T)")
    ax.bar([i+w/2 for i in x], [recu.get(p, 0) for p in patients], w, color=COL_M, label="Recurrent (M)")
    ax.set_xticks(list(x)); ax.set_xticklabels(patients, fontsize=9)
    ax.set_ylabel("PASS variants"); ax.set_title("Total somatic mutation burden")
    ax.legend(); ax.spines[["top","right"]].set_visible(False); fig.tight_layout()
    _save(fig, out_dir / "plot_mutation_burden.png")


# ═════════════════════════════════════════════════════════════════════════════
# 3. FS/SG BURDEN — frameshift + stop_gained per sample
# ═════════════════════════════════════════════════════════════════════════════

def plot_fs_sg_burden(summary: pd.DataFrame, out_dir: Path):
    prim = summary[summary.timepoint == "primary"].set_index("patient")["pass_fs_sg"]
    recu = summary[summary.timepoint == "recurrent"].set_index("patient")["pass_fs_sg"]
    patients = sorted(prim.index.union(recu.index), key=int)
    x, w = range(len(patients)), 0.4
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar([i-w/2 for i in x], [prim.get(p, 0) for p in patients], w, color=COL_T, label="Primary (T)")
    ax.bar([i+w/2 for i in x], [recu.get(p, 0) for p in patients], w, color=COL_M, label="Recurrent (M)")
    ax.set_xticks(list(x)); ax.set_xticklabels(patients, fontsize=9)
    ax.set_ylabel("FS/SG variants (PASS)"); ax.set_title("High-impact truncating variants (FS/SG)")
    ax.legend(); ax.spines[["top","right"]].set_visible(False); fig.tight_layout()
    _save(fig, out_dir / "plot_fs_sg_burden.png")


# ═════════════════════════════════════════════════════════════════════════════
# 4. PAIRED PRIMARY vs RECURRENT — scatter plot per patient
# ═════════════════════════════════════════════════════════════════════════════

def plot_paired_scatter(summary: pd.DataFrame, out_dir: Path):
    prim = summary[summary.timepoint == "primary"].set_index("patient")["pass_fs_sg"]
    recu = summary[summary.timepoint == "recurrent"].set_index("patient")["pass_fs_sg"]
    patients = sorted(prim.index.intersection(recu.index), key=int)
    xv, yv = [prim[p] for p in patients], [recu[p] for p in patients]
    mx = max(max(xv), max(yv)) * 1.1
    above = sum(y > x for x, y in zip(xv, yv))
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot([0, mx], [0, mx], "k--", lw=0.8, alpha=0.4)
    for p, x, y in zip(patients, xv, yv):
        ax.scatter(x, y, color=COL_M if y > x else COL_T, s=60, zorder=3)
        ax.annotate(p, (x, y), xytext=(5, 3), textcoords="offset points", fontsize=8, color="#555")
    ax.set_xlabel("Primary FS/SG"); ax.set_ylabel("Recurrent FS/SG")
    ax.set_title(f"Paired FS/SG burden ({above}/{len(patients)} patients show increase at recurrence)")
    ax.legend(handles=[mpatches.Patch(color=COL_M, label="Higher at recurrence"),
                       mpatches.Patch(color=COL_T, label="Higher at primary")], fontsize=9)
    ax.spines[["top","right"]].set_visible(False); fig.tight_layout()
    _save(fig, out_dir / "plot_paired_scatter.png")


# ═════════════════════════════════════════════════════════════════════════════
# 5. SNV vs INDEL BREAKDOWN — per sample stacked bars
# ═════════════════════════════════════════════════════════════════════════════

def plot_snv_indel_ratio(summary: pd.DataFrame, out_dir: Path):
    patients = sorted(summary["patient"].unique(), key=int)
    t = summary[summary.timepoint == "primary"].set_index("patient")
    m = summary[summary.timepoint == "recurrent"].set_index("patient")
    x, w = range(len(patients)), 0.38
    fig, ax = plt.subplots(figsize=(16, 5))
    for i, p in enumerate(patients):
        for df, sign, c in [(t, -1, COL_T), (m, 1, COL_M)]:
            if p in df.index:
                ax.bar(i+sign*w/2, df.loc[p,"pass_snv"], w, color=c, alpha=0.9)
                ax.bar(i+sign*w/2, df.loc[p,"pass_indel"], w, color=c, alpha=0.45,
                       bottom=df.loc[p,"pass_snv"], hatch="///")
    ax.set_xticks(list(x)); ax.set_xticklabels(patients, fontsize=9)
    ax.set_ylabel("PASS variant count"); ax.set_title("SNV vs indel breakdown")
    ax.legend(handles=[mpatches.Patch(color=COL_T, alpha=0.9, label="Primary SNV"),
                       mpatches.Patch(color=COL_T, alpha=0.45, hatch="///", label="Primary indel"),
                       mpatches.Patch(color=COL_M, alpha=0.9, label="Recurrent SNV"),
                       mpatches.Patch(color=COL_M, alpha=0.45, hatch="///", label="Recurrent indel")], fontsize=9)
    ax.spines[["top","right"]].set_visible(False); fig.tight_layout()
    _save(fig, out_dir / "plot_snv_indel.png")


# ═════════════════════════════════════════════════════════════════════════════
# 6. VARIANT OVERLAP — shared / T-only / M-only per patient
# ═════════════════════════════════════════════════════════════════════════════

def plot_variant_overlap(overlap_df: pd.DataFrame, out_dir: Path):
    patients, x, w = overlap_df["patient"].tolist(), range(len(overlap_df)), 0.55
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar(x, overlap_df["shared"], w, color="#888", label="Shared")
    ax.bar(x, overlap_df["t_only"], w, color=COL_T, bottom=overlap_df["shared"], label="Primary-only")
    ax.bar(x, overlap_df["m_only"], w, color=COL_M,
           bottom=overlap_df["shared"]+overlap_df["t_only"], label="Recurrent-only")
    ax.set_xticks(list(x)); ax.set_xticklabels(patients, fontsize=9)
    ax.set_ylabel("PASS variant count"); ax.set_title("Variant overlap — shared vs timepoint-specific")
    ax.legend(); ax.spines[["top","right"]].set_visible(False); fig.tight_layout()
    _save(fig, out_dir / "plot_variant_overlap.png")


# ═════════════════════════════════════════════════════════════════════════════
# 7. RECURRENCE-ACQUIRED FS/SG — the TMZ neoantigen candidate pool
# ═════════════════════════════════════════════════════════════════════════════

def plot_m_only_fs_sg(overlap_df: pd.DataFrame, out_dir: Path):
    df = overlap_df.sort_values("m_only_fs_sg", ascending=False)
    median = df["m_only_fs_sg"].median()
    colors = [COL_HI if v > median*2 else COL_M for v in df["m_only_fs_sg"]]
    fig, ax = plt.subplots(figsize=(14, 5))
    bars = ax.bar(range(len(df)), df["m_only_fs_sg"], color=colors, width=0.6)
    for bar, pct in zip(bars, df["m_only_fs_sg_pct"]):
        if bar.get_height() > 0:
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
                    f"{pct}%", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(range(len(df))); ax.set_xticklabels(df["patient"].tolist(), fontsize=9)
    ax.set_ylabel("Recurrent-only FS/SG variants")
    ax.set_title("Recurrence-acquired FS/SG per patient (% = fraction of M-only variants)")
    ax.spines[["top","right"]].set_visible(False); fig.tight_layout()
    _save(fig, out_dir / "plot_m_only_fs_sg.png")


# ═════════════════════════════════════════════════════════════════════════════
# 8. GENE-LEVEL RECURRENCE — genes hit by FS/SG across patients
# ═════════════════════════════════════════════════════════════════════════════

def analyse_gene_recurrence(vcf_dir: Path, out_dir: Path) -> pd.DataFrame:
    gene_hits = {}
    for vcf_path in _list_vcfs(vcf_dir):
        m = VCF_PATTERN.match(vcf_path.name)
        if not m:
            continue
        pid, tp = m.group(1), m.group(2)
        for v in VCF(str(vcf_path)):
            if v.FILTER and v.FILTER != "PASS":
                continue
            ann = str(v.INFO.get("ANN", ""))
            if not HIGH_IMPACT.search(ann):
                continue
            gene = ann.split(",")[0].split("|")[3] if ann else "Unknown"
            gene_hits.setdefault(gene, {"T": set(), "M": set()})[tp].add(pid)
    df = pd.DataFrame([{"gene": g, "n_primary": len(h["T"]), "n_recurrent": len(h["M"]),
                        "n_either": len(h["T"]|h["M"]), "recurrent_only": len(h["M"]-h["T"])}
                       for g, h in gene_hits.items()]).sort_values("n_recurrent", ascending=False)
    df.to_csv(out_dir / "gene_recurrence.tsv", sep="\t", index=False)
    return df


def plot_top_recurrent_genes(gene_df: pd.DataFrame, out_dir: Path, top_n: int = 20):
    df = gene_df.head(top_n).sort_values("n_recurrent")
    colors = [COL_HI if (r.recurrent_only > 0 and r.n_primary == 0) else COL_M
              for _, r in df.iterrows()]
    fig, ax = plt.subplots(figsize=(8, top_n*0.45+1.5))
    ax.barh(df["gene"], df["n_recurrent"], color=colors, height=0.6, label="Recurrent (M)")
    ax.barh(df["gene"], df["n_primary"], color=COL_T, height=0.6, alpha=0.6, label="Primary (T)")
    ax.set_xlabel("Patients with FS/SG hit")
    ax.set_title(f"Top {top_n} recurrently mutated genes (amber = recurrent-only)")
    ax.legend(fontsize=9); ax.spines[["top","right"]].set_visible(False); fig.tight_layout()
    _save(fig, out_dir / "plot_top_recurrent_genes.png")


# ═════════════════════════════════════════════════════════════════════════════
# 9. TMZ MUTATIONAL SIGNATURE — C>T at CpG (SBS11 proxy) in recurrent tumors
#    Requires indexed GRCh38 FASTA (.fai). Auto-detected from data/reference/.
# ═════════════════════════════════════════════════════════════════════════════

def _fai_get_context(fai_index: dict, fasta_path: Path, chrom: str, pos: int):
    """Fetch 3-base context around pos (1-based) using .fai byte index."""
    if chrom not in fai_index:
        chrom = chrom.replace("chr","") if chrom.startswith("chr") else "chr"+chrom
    if chrom not in fai_index:
        return None
    info = fai_index[chrom]
    start = pos - 2  # 0-based
    if start < 0:
        return None
    lines, rem = divmod(start, info["bases"])
    offset = info["offset"] + lines * info["bytes"] + rem
    with open(fasta_path, "rb") as f:
        f.seek(offset)
        raw = f.read(12).decode("ascii", errors="ignore")
    bases = "".join(c for c in raw if c.isalpha())[:3]
    return bases.upper() if len(bases) == 3 else None


def compute_tmz_signature(vcf_dir: Path, fasta_path: Path) -> pd.DataFrame:
    fai_path = Path(str(fasta_path) + ".fai")
    if not fasta_path.exists() or not fai_path.exists():
        print(f"[WARN] TMZ signature: FASTA/index not found at {fasta_path}")
        return pd.DataFrame()
    # Load .fai index
    fai = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                fai[parts[0]] = {"offset": int(parts[2]), "bases": int(parts[3]), "bytes": int(parts[4])}
    complement = {"A":"T","T":"A","C":"G","G":"C"}
    records = []
    for vcf_path in _list_vcfs(vcf_dir):
        m = VCF_PATTERN.match(vcf_path.name)
        if not m or m.group(2) != "M":
            continue
        pid = int(m.group(1))
        total = tmz = 0
        for v in VCF(str(vcf_path)):
            if v.FILTER and v.FILTER != "PASS":
                continue
            if len(v.REF) != 1 or len(v.ALT[0]) != 1:
                continue
            total += 1
            ref, alt = v.REF.upper(), v.ALT[0].upper()
            r, a = (ref, alt) if ref in ("C","T") else (complement.get(ref,ref), complement.get(alt,alt))
            if r == "C" and a == "T":
                ctx = _fai_get_context(fai, fasta_path, v.CHROM, v.POS)
                if ctx and ctx[2] == "G":
                    tmz += 1
        pct = round(tmz/total*100, 1) if total else 0.0
        records.append({"patient": pid, "total_snv": total, "tmz_snv": tmz, "tmz_pct": pct})
        print(f"  Patient {pid}: {tmz}/{total} C>T@CpG ({pct}%)")
    return pd.DataFrame(records).sort_values("patient")


def plot_tmz_signature(tmz_df: pd.DataFrame, out_dir: Path):
    df = tmz_df.sort_values("tmz_pct", ascending=False)
    median = df["tmz_pct"].median()
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar(range(len(df)), df["tmz_pct"],
           color=[COL_M if v >= median else "#e0c070" for v in df["tmz_pct"]], width=0.6)
    ax.axhline(median, color="#555", linestyle="--", lw=1, label=f"Median: {median:.1f}%")
    ax.set_xticks(range(len(df))); ax.set_xticklabels([str(int(p)) for p in df["patient"]], fontsize=9)
    ax.set_ylabel("C>T at CpG (% of PASS SNVs)")
    ax.set_title("TMZ mutational signature — recurrent tumors (SBS11 proxy)")
    ax.legend(); ax.spines[["top","right"]].set_visible(False); fig.tight_layout()
    _save(fig, out_dir / "plot_tmz_signature.png")


# ═════════════════════════════════════════════════════════════════════════════
# 10. NMD CANDIDATE PRIORITISATION — ranking patients by FS/SG candidate pool
#     Computed inside generate_html_report from overlap_df
# ═════════════════════════════════════════════════════════════════════════════


# ═════════════════════════════════════════════════════════════════════════════
# HLA TYPING (Stage 4) — load Dragen HLA TSVs
# ═════════════════════════════════════════════════════════════════════════════

HLA_GENES = ["A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1"]

def load_hla_typing(hla_dir: Path) -> pd.DataFrame:
    """Load HLA Class I alleles from CSV (sample_id, HLA_Types format)."""
    csv = list(hla_dir.glob("*.csv"))
    if csv:
        df = pd.read_csv(csv[0])
        print(f"[INFO] HLA CSV loaded: {len(df)} samples from {csv[0].name}")
        return df
    # Fallback: Dragen TSV files
    tsv_pat = re.compile(r"^S(\d+)_(T|M)\.hla\.tsv$")
    records = []
    for f in sorted(hla_dir.iterdir()):
        m = tsv_pat.match(f.name)
        if not m:
            continue
        sample = f"{m.group(1)}{m.group(2)}"
        df = pd.read_csv(f, sep="\t")
        alleles = []
        for gene in ["A","B","C"]:
            row = df[df["gene"] == gene]
            if row.empty:
                continue
            for col in ["allele_1","allele_2"]:
                a = row.iloc[0].get(col)
                if pd.notna(a) and str(a) not in ("NA",""):
                    alleles.append(f"HLA-{gene}*{a}" if "*" not in str(a) else f"HLA-{a}")
        records.append({"Sample_ID": sample, "HLA_Types": ";".join(alleles)})
    return pd.DataFrame(records)


# ═════════════════════════════════════════════════════════════════════════════
# HTML REPORT GENERATOR
# ═════════════════════════════════════════════════════════════════════════════

def generate_html_report(out_dir, summary, overlap_df, gene_df, hla_df,
                         tmz_df, run_ts) -> Path:
    import base64

    def b64(name):
        p = out_dir / name
        return ("data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()) if p.exists() else ""

    def fig(src, caption=""):
        if not src: return ""
        return (f'<div style="margin:16px 0;"><img src="{src}" style="max-width:75%;display:block;'
                f'margin:0 auto;border-radius:4px;box-shadow:0 1px 4px rgba(0,0,0,.12);">'
                + (f'<p style="font-size:11px;color:#888;margin-top:6px;">{caption}</p>' if caption else "")
                + "</div>")

    def sec(title, body):
        return (f'<section style="margin-bottom:40px;"><h2 style="font-size:16px;font-weight:600;'
                f'color:#1a1a1a;border-bottom:1px solid #e5e5e5;padding-bottom:6px;margin-bottom:16px;">'
                f'{title}</h2>{body}</section>')

    def card(label, value, sub=""):
        return (f'<div style="background:#f8f8f8;border:1px solid #e5e5e5;border-radius:6px;'
                f'padding:12px 16px;min-width:120px;"><div style="font-size:11px;color:#888;">{label}</div>'
                f'<div style="font-size:20px;font-weight:600;">{value}</div>'
                + (f'<div style="font-size:10px;color:#aaa;margin-top:2px;">{sub}</div>' if sub else "")
                + "</div>")

    def cards(*items):
        return '<div style="display:flex;gap:10px;flex-wrap:wrap;margin-bottom:20px;">' + "".join(items) + "</div>"

    def tbl(headers, rows):
        th = "".join(f'<th style="background:#f4f4f4;padding:5px 10px;text-align:{"left" if i==0 else "right"};'
                     f'font-weight:500;border-bottom:2px solid #ddd;">{h}</th>'
                     for i, h in enumerate(headers))
        tr = "".join("<tr>" + "".join(
            f'<td style="padding:4px 10px;border-bottom:1px solid #eee;text-align:{"left" if i==0 else "right"};'
            f'font-size:12px;">{v}</td>' for i, v in enumerate(row)) + "</tr>" for row in rows)
        return (f'<div style="overflow-x:auto;margin-top:10px;"><table style="width:100%;border-collapse:collapse;">'
                f'<thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table></div>')

    # ── compute stats ─────────────────────────────────────────────────────────
    n_p   = int((summary.timepoint=="primary").sum())
    n_m   = int((summary.timepoint=="recurrent").sum())
    fs_t  = int(summary[summary.timepoint=="primary"]["pass_fs_sg"].sum())
    fs_m  = int(summary[summary.timepoint=="recurrent"]["pass_fs_sg"].sum())
    pct_i = round((fs_m-fs_t)/fs_t*100,1) if fs_t else 0
    med_t = round(float(summary[summary.timepoint=="primary"]["pass_fs_sg"].median()),1)
    med_m = round(float(summary[summary.timepoint=="recurrent"]["pass_fs_sg"].median()),1)
    top_s = summary.loc[summary.pass_fs_sg.idxmax(),"sample"]
    top_n = int(summary.pass_fs_sg.max())
    hm    = sorted(summary[summary.pass_total>HYPERMUTATOR_THRESHOLD]["sample"].tolist())

    # section 10: NMD prioritisation
    prio = overlap_df.sort_values("m_only_fs_sg",ascending=False).reset_index(drop=True)
    prio["prio_rank"] = range(1,len(prio)+1)
    q75, q50 = prio["m_only_fs_sg"].quantile(0.75), prio["m_only_fs_sg"].quantile(0.50)
    prio["priority"] = prio["m_only_fs_sg"].apply(
        lambda n: "High" if n>=q75 else ("Medium" if n>=q50 else "Low"))

    # TMZ stats
    if tmz_df is not None and not tmz_df.empty:
        tmz_med   = round(float(tmz_df["tmz_pct"].median()),1)
        tmz_tot   = int(tmz_df["total_snv"].sum())
        tmz_hits  = int(tmz_df["tmz_snv"].sum())
        tmz_pct   = round(tmz_hits/tmz_tot*100,1) if tmz_tot else 0
    else:
        tmz_med = tmz_tot = tmz_hits = tmz_pct = 0

    hla_note = f"{len(hla_df)} sample(s) loaded" if hla_df is not None else "Not provided"
    top5     = overlap_df.nlargest(5,"m_only_fs_sg")
    top_genes= gene_df.head(5)

    # ── load plot images ──────────────────────────────────────────────────────
    B = {k: b64(v) for k, v in {
        "burden":"plot_mutation_burden.png", "fs":"plot_fs_sg_burden.png",
        "paired":"plot_paired_scatter.png", "snv":"plot_snv_indel.png",
        "overlap":"plot_variant_overlap.png", "monly":"plot_m_only_fs_sg.png",
        "genes":"plot_top_recurrent_genes.png", "tmz":"plot_tmz_signature.png"}.items()}

    # ── helper for per-patient lookup ─────────────────────────────────────────
    def sv(patient, tp, col):
        r = summary[(summary.patient==patient)&(summary.timepoint==tp)]
        return f"{int(r[col].values[0]):,}" if len(r) else "-"

    CSS = ("* {box-sizing:border-box;margin:0;padding:0;}"
           "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;"
           "font-size:13px;color:#333;background:#fafafa;}"
           ".wrap{max-width:960px;margin:0 auto;padding:32px 24px;}"
           "header{margin-bottom:40px;padding-bottom:20px;border-bottom:2px solid #e5e5e5;}"
           "header h1{font-size:22px;font-weight:600;} header p{font-size:12px;color:#888;margin-top:4px;}"
           "img{max-width:100%;height:auto;}")

    patients = sorted(summary.patient.unique())

    sections = [
        # ── 1. Dataset overview ───────────────────────────────────────────────
        sec("1. Dataset overview",
            f"<p style='margin-bottom:14px;'>Somatic variant calls from WES. {n_p} primary (T) and "
            f"{n_m} recurrent (M) samples. PASS-filtered variants only.</p>"
            + cards(card("Samples",len(summary),f"{min(n_p,n_m)} matched pairs"),
                    card("Primary (T)",n_p,"pre-treatment"),
                    card("Recurrent (M)",n_m,"post-TMZ"),
                    card("HLA typed",len(hla_df) if hla_df is not None else 0,hla_note))),

        # ── 2. Total mutation burden ──────────────────────────────────────────
        sec("2. Total somatic mutation burden",
            ("<p style='margin-bottom:12px;'>Total PASS variant counts per sample."
             + (f" Potential hypermutators: <strong>{', '.join(hm)}</strong> (>{HYPERMUTATOR_THRESHOLD} variants)." if hm else "")
             + "</p>")
            + fig(B["burden"],"Fig 1. Total somatic mutation burden.")
            + tbl(["Patient","Primary PASS","Recurrent PASS","Δ"],
                  [[p, sv(p,"primary","pass_total"), sv(p,"recurrent","pass_total"),
                    f"+{int(summary[(summary.patient==p)&(summary.timepoint=='recurrent')]['pass_total'].values[0])-int(summary[(summary.patient==p)&(summary.timepoint=='primary')]['pass_total'].values[0]):,}"
                    if len(summary[(summary.patient==p)&(summary.timepoint=='primary')])>0 and len(summary[(summary.patient==p)&(summary.timepoint=='recurrent')])>0 else "-"]
                   for p in patients])),

        # ── 3. FS/SG burden ───────────────────────────────────────────────────
        sec("3. High-impact truncating variants (FS/SG)",
            f"<p style='margin-bottom:12px;'>Frameshift and stop-gained variants — the mutation class most "
            f"relevant to NMD. Recurrent tumors show a {pct_i}% increase (median: {med_t} → {med_m}/sample). "
            f"Top outlier: {top_s} ({top_n:,} variants).</p>"
            + fig(B["fs"],"Fig 2. FS/SG burden per sample.")
            + tbl(["Patient","Primary FS/SG","Recurrent FS/SG","Δ"],
                  [[p, sv(p,"primary","pass_fs_sg"), sv(p,"recurrent","pass_fs_sg"),
                    f"+{int(summary[(summary.patient==p)&(summary.timepoint=='recurrent')]['pass_fs_sg'].values[0])-int(summary[(summary.patient==p)&(summary.timepoint=='primary')]['pass_fs_sg'].values[0]):,}"
                    if len(summary[(summary.patient==p)&(summary.timepoint=='primary')])>0 and len(summary[(summary.patient==p)&(summary.timepoint=='recurrent')])>0 else "-"]
                   for p in patients])),

        # ── 4. Paired comparison ──────────────────────────────────────────────
        sec("4. Paired primary vs recurrent comparison",
            "<p style='margin-bottom:12px;'>Direct paired comparison of FS/SG burden per patient.</p>"
            + fig(B["paired"],"Fig 3. Paired scatter: primary vs recurrent FS/SG burden.")
            + tbl(["Patient","Primary","Recurrent","M-only FS/SG","% FS/SG"],
                  [[str(int(r.patient)),f"{int(r.t_total):,}",f"{int(r.m_total):,}",
                    f"{int(r.m_only_fs_sg):,}",f"{r.m_only_fs_sg_pct:.1f}%"]
                   for _,r in overlap_df.iterrows()])),

        # ── 5. SNV vs indel ───────────────────────────────────────────────────
        sec("5. SNV vs indel breakdown",
            "<p style='margin-bottom:12px;'>Indels include frameshift insertions/deletions — "
            "the NMD-relevant class. The TMZ signature is predominantly C→T SNVs.</p>"
            + fig(B["snv"],"Fig 4. SNV vs indel per sample.")
            + tbl(["Patient","Timepoint","SNV","Indel","% Indel"],
                  [[r.patient,r.timepoint,f"{int(r.pass_snv):,}",f"{int(r.pass_indel):,}",
                    f"{round(r.pass_indel/(r.pass_snv+r.pass_indel)*100,1)}%"
                    if (r.pass_snv+r.pass_indel)>0 else "-"]
                   for _,r in summary.iterrows()])),

        # ── 6. Variant overlap ────────────────────────────────────────────────
        sec("6. Variant overlap — shared vs timepoint-specific",
            "<p style='margin-bottom:12px;'>Variants shared between primary and recurrent vs unique to each.</p>"
            + fig(B["overlap"],"Fig 5. Variant overlap per patient.")
            + tbl(["Patient","Shared","Primary-only","Recurrent-only","M-only FS/SG"],
                  [[str(int(r.patient)),f"{int(r.shared):,}",f"{int(r.t_only):,}",
                    f"{int(r.m_only):,}",f"{int(r.m_only_fs_sg):,}"]
                   for _,r in overlap_df.iterrows()])),

        # ── 7. Recurrence-acquired FS/SG ──────────────────────────────────────
        sec("7. Recurrence-acquired FS/SG — TMZ candidate pool",
            "<p style='margin-bottom:12px;'>FS/SG variants unique to the recurrent tumor — "
            "the primary NMD-neoantigen candidates.</p>"
            + fig(B["monly"],"Fig 6. Recurrence-acquired FS/SG per patient.")
            + tbl(["Patient","Recurrent-only","FS/SG","% FS/SG"],
                  [[str(int(r.patient)),f"{int(r.m_only):,}",
                    f"{int(r.m_only_fs_sg):,}",f"{r.m_only_fs_sg_pct:.1f}%"]
                   for _,r in top5.iterrows()])
            + "<p style='font-size:11px;color:#888;margin-top:8px;'>Top 5 shown. Full data in paired_variant_overlap.tsv.</p>"),

        # ── 8. Gene-level recurrence ──────────────────────────────────────────
        sec("8. Gene-level recurrence",
            "<p style='margin-bottom:12px;'>Genes recurrently hit by FS/SG across patients. "
            "Recurrent-only genes are TMZ-induced neoantigen candidates.</p>"
            + fig(B["genes"],"Fig 7. Top recurrently mutated genes.")
            + tbl(["Gene","Primary patients","Recurrent patients","Recurrent-only"],
                  [[r.gene,r.n_primary,r.n_recurrent,r.recurrent_only]
                   for _,r in top_genes.iterrows()])),

        # ── 9. TMZ signature ──────────────────────────────────────────────────
        sec("9. TMZ mutational signature — SBS11 enrichment in recurrent tumors",
            "<p style='margin-bottom:12px;'>TMZ causes C>T transitions at CpG dinucleotides (COSMIC SBS11). "
            "Enrichment in recurrent tumors confirms variants are therapy-induced.</p>"
            + (cards(card("Cohort C>T@CpG",f"{tmz_pct}%","of recurrent SNVs"),
                     card("Median per patient",f"{tmz_med}%","C>T at CpG"))
               + fig(B["tmz"],"Fig 8. C>T at CpG per patient (recurrent tumors).")
               + f"<p style='font-size:12px;color:#555;'>{tmz_hits:,} of {tmz_tot:,} recurrent PASS SNVs "
               f"({tmz_pct}%) match the TMZ C>T@CpG signature.</p>"
               if tmz_tot > 0
               else "<p style='color:#888;'>TMZ analysis unavailable — ensure GRCh38 FASTA is indexed in data/reference/.</p>")),

        # ── 10. NMD candidate prioritisation ─────────────────────────────────
        sec("10. NMD neoantigen candidate prioritisation",
            "<p style='margin-bottom:12px;'>Patients ranked by recurrence-acquired FS/SG variants — "
            "the pool from which NMD-silenced neoantigens are drawn. "
            "<strong>High</strong> tier (top 25%) are primary candidates for pVACseq and NMD analysis.</p>"
            + tbl(["Rank","Patient","Recurrence-acquired FS/SG","% of M-only","Priority"],
                  [[f"#{int(r.prio_rank)}", str(int(r.patient)), str(int(r.m_only_fs_sg)),
                    f"{r.m_only_fs_sg_pct:.1f}%",
                    (f'<span style="background:{"#e8f5e9" if r.priority=="High" else "#fff8e1" if r.priority=="Medium" else "#f5f5f5"};'
                     f'color:{"#2e7d32" if r.priority=="High" else "#e65100" if r.priority=="Medium" else "#616161"};'
                     f'padding:2px 8px;border-radius:10px;font-size:11px;font-weight:600;">{r.priority}</span>')]
                   for _,r in prio.iterrows()])
            + f"<p style='font-size:11px;color:#888;margin-top:8px;'>"
            f"High ≥{q75:.0f} FS/SG; Medium ≥{q50:.0f}; Low = bottom 50%. Full data in paired_variant_overlap.tsv.</p>"),
    ]

    html = (f'<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">'
            f'<title>GBM Mutation Landscape — {run_ts}</title><style>{CSS}</style></head>'
            f'<body><div class="wrap"><header><h1>GBM Somatic Mutation Landscape</h1>'
            f'<p>Run: {run_ts} &middot; GBM NMD-Neoantigen Pipeline &middot; '
            f'github.com/paleslui/gbm-nmd-pipeline</p></header>'
            + "\n".join(sections) + "</div></body></html>")

    out = out_dir / "report.html"
    out.write_text(html, encoding="utf-8")
    print(f"[REPORT] Saved {out}")
    return out


# ═════════════════════════════════════════════════════════════════════════════
# UTILITY
# ═════════════════════════════════════════════════════════════════════════════

def _save(fig, path: Path):
    fig.savefig(path, dpi=150); plt.close(fig)
    print(f"[PLOT] Saved {path}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="GBM NMD-Neoantigen Pipeline — Stages 1-4")
    parser.add_argument("--vcf_dir",  required=True, help="Directory of SnpEff-annotated VCFs")
    parser.add_argument("--out_dir",  default="./gbm_output", help="Output directory")
    parser.add_argument("--hla_dir",  default=None, help="HLA typing directory (Stage 4)")
    parser.add_argument("--fasta",    default=None, help="GRCh38 FASTA for TMZ signature (optional)")
    args = parser.parse_args()

    vcf_dir  = Path(args.vcf_dir)
    # Output directly to --out_dir (master_pipeline.sh provides run_<TS> at the parent level)
    run_ts   = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir  = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta    = Path(args.fasta) if args.fasta else vcf_dir.parent / "reference" / "GRCh38.primary_assembly.fa"

    print(f"\n{'='*60}\n  GBM NMD-Neoantigen Pipeline — Stages 1-4\n{'='*60}\n")

    # Stage 1-2: parse VCFs
    print("[STAGE 1-2] Parsing VCFs...")
    summary, fs_sg_dict = build_summary(vcf_dir)
    summary.to_csv(out_dir/"summary_mutation_burden.tsv", sep="\t", index=False)
    pd.concat(fs_sg_dict.values()).to_csv(out_dir/"all_fs_sg_variants.tsv", sep="\t", index=False)

    # Stage 3: plots
    print("\n[STAGE 3] Generating plots...")
    plot_mutation_burden(summary, out_dir)
    plot_fs_sg_burden(summary, out_dir)
    plot_paired_scatter(summary, out_dir)
    plot_snv_indel_ratio(summary, out_dir)

    # Stage 3B: paired variant analysis
    print("\n[STAGE 3B] Comparing paired variants...")
    overlap_df = compare_paired_variants(vcf_dir)
    overlap_df.to_csv(out_dir/"paired_variant_overlap.tsv", sep="\t", index=False)
    plot_variant_overlap(overlap_df, out_dir)
    plot_m_only_fs_sg(overlap_df, out_dir)

    # Stage 3C: gene recurrence
    print("\n[STAGE 3C] Gene-level recurrence...")
    gene_df = analyse_gene_recurrence(vcf_dir, out_dir)
    plot_top_recurrent_genes(gene_df, out_dir)

    # Stage 3D: TMZ signature
    print("\n[STAGE 3D] TMZ mutational signature...")
    tmz_df = compute_tmz_signature(vcf_dir, fasta)
    if not tmz_df.empty:
        tmz_df.to_csv(out_dir/"tmz_signature.tsv", sep="\t", index=False)
        plot_tmz_signature(tmz_df, out_dir)

    # Stage 4: HLA typing
    hla_df = None
    if args.hla_dir:
        print("\n[STAGE 4] Loading HLA typing...")
        try:
            hla_df = load_hla_typing(Path(args.hla_dir))
            hla_df.to_csv(out_dir/"hla_typing_summary.tsv", sep="\t", index=False)
        except Exception as e:
            print(f"  [WARN] HLA loading failed: {e}")
    else:
        print("\n[STAGE 4] HLA typing skipped (no --hla_dir provided)")

    # HTML report
    print("\n[REPORT] Generating HTML report...")
    generate_html_report(out_dir, summary, overlap_df, gene_df, hla_df, tmz_df, run_ts)

    print(f"\n{'='*60}\n  Done. Output: {out_dir}\n{'='*60}\n")


if __name__ == "__main__":
    main()