#!/usr/bin/env python3
"""
nmd_scoring.py — GBM Pipeline Stage 7: NMD Sensitivity Scoring
Part of the GBM NMD-Neoantigen Pipeline.
https://github.com/paleslui/gbm-nmd-pipeline

Scores each pVACseq neoantigen candidate for NMD sensitivity using:
  1. VEP NMD plugin (CSQ field) — requires Dragen VEP-annotated VCF
  2. Lindeboom et al. 2019 rule-based method (4 rules)

Ensemble confidence 0-3: 3=both agree, 2=single method, 1=disagree, 0=no data

Usage:
  python nmd_scoring.py --pvacseq_tsv <path> --vep_vcf <path> --out_dir <path>
  python nmd_scoring.py --pvacseq_tsv <path> --out_dir <path>  # rule-based only
"""

import argparse, gzip, re
from pathlib import Path

import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io, base64

COL_S = "#D85A30"   # sensitive  — coral
COL_I = "#378ADD"   # insensitive — blue
COL_U = "#F0C040"   # uncertain  — amber

TRUNCATING = ['stop_gained', 'frameshift', 'stop_lost', 'start_lost',
              'splice_donor', 'splice_acceptor']


# ═════════════════════════════════════════════════════════════════════════════
# VEP VCF PARSING
#   parse_csq_header: extract CSQ field names from VCF header
#   parse_vep_vcf: return dict of (chrom,pos,ref,alt) → transcript list
# ═════════════════════════════════════════════════════════════════════════════

def parse_csq_header(vcf_path: Path) -> list:
    opener = gzip.open if str(vcf_path).endswith('.gz') else open
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('##INFO=<ID=CSQ'):
                m = re.search(r'Format: ([^"]+)"', line)
                return m.group(1).strip().split('|') if m else []
    return []


def parse_vep_vcf(vcf_path: Path) -> dict:
    fields = parse_csq_header(vcf_path)
    if not fields:
        print("[WARN] Could not parse CSQ header — VEP NMD plugin scoring skipped")
        return {}
    print(f"[INFO] CSQ fields found: {len(fields)}")
    print(f"[INFO] NMD plugin field present: {'NMD' in fields}")
    print(f"[INFO] CANONICAL field present:  {'CANONICAL' in fields}")
    variants = {}
    opener = gzip.open if str(vcf_path).endswith('.gz') else open
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            chrom, pos, ref, alt = parts[0], int(parts[1]), parts[3], parts[4]
            csq = next((f[4:] for f in parts[7].split(';') if f.startswith('CSQ=')), '')
            if not csq:
                continue
            transcripts = [dict(zip(fields, e.split('|')))
                           for e in csq.split(',') if len(e.split('|')) == len(fields)]
            variants[(chrom, pos, ref, alt)] = transcripts
    print(f"[INFO] Parsed {len(variants)} variants from VEP VCF")
    return variants


# ═════════════════════════════════════════════════════════════════════════════
# METHOD 1: VEP NMD PLUGIN
#   Reads the NMD field from the CSQ annotation.
#   Empty NMD field on truncating variant = SENSITIVE.
#   'NMD_escaping_variant' = INSENSITIVE.
# ═════════════════════════════════════════════════════════════════════════════

def score_nmd_vep_plugin(transcripts: list) -> str:
    if not transcripts:
        return 'UNKNOWN'
    canon = [t for t in transcripts if t.get('CANONICAL') == 'YES'
             and t.get('BIOTYPE') == 'protein_coding']
    t = canon[0] if canon else transcripts[0]
    if 'NMD' not in t:
        return 'UNKNOWN'
    if not any(c in t.get('Consequence', '') for c in TRUNCATING):
        return 'UNKNOWN'
    nmd = t.get('NMD', '').strip()
    if not nmd:
        return 'SENSITIVE'
    return 'INSENSITIVE' if 'escaping' in nmd.lower() else 'UNKNOWN'


# ═════════════════════════════════════════════════════════════════════════════
# METHOD 2: LINDEBOOM RULE-BASED SCORING (Nat Cell Biol 2019)
#   Rules applied in priority order (first match wins):
#     Rule 4: PTC within 150nt of start codon       → INSENSITIVE (pioneer round)
#     Rule 1: PTC in last exon                       → INSENSITIVE (no downstream EJC)
#     Rule 3: PTC in long exon (>407nt)              → INSENSITIVE (EJC too far)
#     Rule 2: PTC >55nt upstream of last EJC         → SENSITIVE (canonical NMD)
# ═════════════════════════════════════════════════════════════════════════════

def _parse_pos(s: str) -> tuple:
    """Parse 'start-end/total' or 'pos/total' → (start, end, total)."""
    if not s:
        return None, None, None
    try:
        main, *rest = s.split('/')
        total = int(rest[0]) if rest else None
        if '-' in main:
            a, b = main.split('-')
            return int(a), int(b), total
        v = int(main)
        return v, v, total
    except (ValueError, IndexError):
        return None, None, None


def _parse_exon(s: str) -> tuple:
    """Parse 'current/total' → (current, total)."""
    if not s or '/' not in s:
        return None, None
    try:
        a, b = s.split('/')
        return int(a), int(b)
    except ValueError:
        return None, None


def score_nmd_rules(transcripts: list) -> tuple:
    """Returns (classification, rule_applied, notes)."""
    if not transcripts:
        return 'UNKNOWN', 'no_transcripts', ''
    canon = [t for t in transcripts if t.get('CANONICAL') == 'YES'
             and t.get('BIOTYPE') == 'protein_coding']
    t = canon[0] if canon else transcripts[0]
    if not any(c in t.get('Consequence', '') for c in TRUNCATING):
        return 'NOT_APPLICABLE', 'not_truncating', ''

    exon_curr, exon_tot   = _parse_exon(t.get('EXON', ''))
    cds_start, _, _       = _parse_pos(t.get('CDS_position', ''))
    cdna_start, cdna_end, _ = _parse_pos(t.get('cDNA_position', ''))

    notes = []
    if exon_curr:  notes.append(f"exon {exon_curr}/{exon_tot}")
    if cds_start:  notes.append(f"CDS {cds_start}nt")
    n = '; '.join(notes)

    # Rule 4: start-proximal PTC
    if cds_start is not None and cds_start <= 150:
        return 'INSENSITIVE', 'rule4_start_proximal_<150nt', n

    if exon_curr is not None and exon_tot is not None:
        # Rule 1: last exon
        if exon_curr == exon_tot:
            return 'INSENSITIVE', 'rule1_last_exon', n
        # Rule 3: long exon >407nt
        if cdna_start is not None and cdna_end is not None:
            exon_len = abs(cdna_end - cdna_start) + 1
            if exon_len > 407:
                return 'INSENSITIVE', f'rule3_long_exon_>407nt', f"{n}; exon ~{exon_len}nt"
        # Rule 2: 55nt boundary (canonical NMD trigger)
        if exon_tot - exon_curr >= 1:
            return 'SENSITIVE', 'rule2_55nt_boundary', n

    return 'UNCERTAIN', 'insufficient_info', n


# ═════════════════════════════════════════════════════════════════════════════
# ENSEMBLE CONSENSUS
#   Combines VEP plugin and rule-based results.
#   Confidence 0-3: 3=both agree, 2=single method, 1=disagree, 0=no data
# ═════════════════════════════════════════════════════════════════════════════

def ensemble_nmd(vep: str, rule: str) -> tuple:
    """Returns (consensus, confidence_label, confidence_score)."""
    v = vep  if vep  in ('SENSITIVE', 'INSENSITIVE') else None
    r = rule if rule in ('SENSITIVE', 'INSENSITIVE') else None
    if v is None and r is None: return 'UNKNOWN',    'no_data',           0
    if v is None:               return r,            'single_method_rules', 2
    if r is None:               return v,            'single_method_vep',   2
    if v == r:                  return v,            'high_confidence',     3
    return                             'UNCERTAIN',  'methods_disagree',    1


# ═════════════════════════════════════════════════════════════════════════════
# MAIN SCORING
#   score_candidates: annotate each pVACseq row with NMD classification
#   hla_allele_breakdown: summarise Tier 1/2 by HLA allele
# ═════════════════════════════════════════════════════════════════════════════

RULE_EXPLANATIONS = {
    'rule1_last_exon':             'PTC in last exon — no downstream EJC',
    'rule2_55nt_boundary':         'PTC >55nt upstream of last exon junction — NMD-sensitive',
    'rule3_long_exon_>407nt':      'PTC in long exon (>407nt) — reduced EJC density',
    'rule4_start_proximal_<150nt': 'PTC within 150nt of start codon — ribosome reinitiation',
}


def score_candidates(pvacseq_df: pd.DataFrame, vep_variants: dict) -> pd.DataFrame:
    results = []
    for _, row in pvacseq_df.iterrows():
        chrom, pos, ref, alt = row['Chromosome'], int(row['Start']), row['Reference'], row['Variant']
        transcripts = vep_variants.get((chrom, pos, ref, alt),
                      vep_variants.get((chrom, pos+1, ref, alt), []))

        vep_nmd                       = score_nmd_vep_plugin(transcripts)
        rule_nmd, rule_applied, notes = score_nmd_rules(transcripts)
        consensus, confidence, score  = ensemble_nmd(vep_nmd, rule_nmd)

        best_ic50   = float(row.get('Best MT IC50 Score', 9999))
        median_ic50 = float(row.get('Median MT IC50 Score', 9999))

        tier = ('TIER1'         if consensus == 'SENSITIVE'   and best_ic50 < 50  else
                'TIER2'         if consensus == 'SENSITIVE'   and best_ic50 < 500 else
                'TIER3_control' if consensus == 'INSENSITIVE' and best_ic50 < 500 else
                'UNCLASSIFIED')

        results.append({**row.to_dict(),
            'nmd_vep_plugin':       vep_nmd,
            'nmd_rules':            rule_nmd,
            'nmd_rule_applied':     rule_applied,
            'nmd_rule_explanation': RULE_EXPLANATIONS.get(rule_applied, rule_applied or 'no transcript data'),
            'nmd_rule_notes':       notes,
            'nmd_consensus':        consensus,
            'nmd_confidence':       confidence,
            'nmd_confidence_score': score,
            'hla_allele':           str(row.get('HLA Allele', '')),
            'best_ic50':            best_ic50,
            'median_ic50':          median_ic50,
            'best_ic50_method':     str(row.get('Best MT IC50 Score Method', '')),
            'priority_tier':        tier,
        })
    return pd.DataFrame(results)


def hla_allele_breakdown(scored_df: pd.DataFrame) -> pd.DataFrame:
    """Summarise Tier 1 and Tier 2 candidates by HLA allele."""
    t12 = scored_df[scored_df['priority_tier'].isin(['TIER1','TIER2'])]
    if t12.empty:
        return pd.DataFrame()
    return (t12.groupby('hla_allele')
            .agg(n_candidates=('priority_tier','count'),
                 n_tier1=('priority_tier', lambda x: (x=='TIER1').sum()),
                 best_ic50=('best_ic50','min'),
                 median_ic50=('median_ic50','median'))
            .reset_index().sort_values('n_tier1', ascending=False))



# ═════════════════════════════════════════════════════════════════════════════
# VISUALISATIONS
#   _b64fig:        convert matplotlib figure to base64 PNG for HTML embedding
#   plot_tiers:     horizontal bar — candidates per tier
#   plot_ic50:      scatter — IC50 by candidate, colored by NMD class
#   plot_nmd_breakdown: grouped bar — classification per method
#   plot_confidence: bar — confidence score distribution
# ═════════════════════════════════════════════════════════════════════════════

def _b64fig(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def plot_tiers(scored_df: pd.DataFrame) -> str:
    """Horizontal bar: candidate counts per priority tier."""
    tiers  = ['TIER1', 'TIER2', 'TIER3_control', 'UNCLASSIFIED']
    labels = ['Tier 1\nNMD-sensitive\nIC50<50nM', 'Tier 2\nNMD-sensitive\nIC50<500nM',
              'Tier 3\nNMD-insensitive\ncontrol', 'Unclassified']
    colors = [COL_S, '#e87040', COL_I, '#cccccc']
    counts = [int((scored_df['priority_tier'] == t).sum()) for t in tiers]

    fig, ax = plt.subplots(figsize=(8, 3))
    bars = ax.barh(labels[::-1], counts[::-1], color=colors[::-1], height=0.5)
    for bar, n in zip(bars, counts[::-1]):
        if n > 0:
            ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2,
                    str(n), va='center', fontsize=11, fontweight='bold')
    ax.set_xlabel("Number of candidates")
    ax.set_title("Candidate count by priority tier", fontsize=12)
    ax.spines[["top","right","left"]].set_visible(False)
    ax.set_xlim(0, max(counts) * 1.3 + 1)
    fig.tight_layout()
    return _b64fig(fig)


def plot_ic50(scored_df: pd.DataFrame) -> str:
    """Scatter: IC50 per candidate colored by NMD consensus."""
    color_map = {'SENSITIVE': COL_S, 'INSENSITIVE': COL_I,
                 'UNCERTAIN': COL_U, 'UNKNOWN': '#cccccc', 'NOT_APPLICABLE': '#aaaaaa'}
    df = scored_df.copy()
    df = df[df['best_ic50'] < 9000].sort_values('best_ic50')
    df['label'] = df['Gene Name'] + '\n' + df['hla_allele'].str.replace('HLA-','')

    fig, ax = plt.subplots(figsize=(10, max(4, len(df)*0.35 + 1)))
    colors = [color_map.get(c, '#aaa') for c in df['nmd_consensus']]
    ax.scatter(df['best_ic50'], range(len(df)), c=colors, s=80, zorder=3)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['label'], fontsize=8)
    ax.axvline(50,  color='#555', linestyle='--', lw=0.8, label='50 nM (Tier 1 threshold)')
    ax.axvline(500, color='#999', linestyle=':', lw=0.8, label='500 nM (Tier 2/3 threshold)')
    ax.set_xlabel("Best MT IC50 Score (nM)")
    ax.set_title("IC50 per candidate — colored by NMD consensus", fontsize=12)
    from matplotlib.patches import Patch
    legend_handles = [Patch(color=v, label=k) for k, v in color_map.items()
                      if k in df['nmd_consensus'].values]
    ax.legend(handles=legend_handles, fontsize=8, loc='lower right')
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    return _b64fig(fig)


def plot_nmd_breakdown(scored_df: pd.DataFrame) -> str:
    """Grouped bar: NMD classification per method (VEP, Rules, Consensus)."""
    cats   = ['SENSITIVE', 'INSENSITIVE', 'NOT_APPLICABLE', 'UNCERTAIN', 'UNKNOWN']
    colors = [COL_S, COL_I, '#aaaaaa', COL_U, '#dddddd']
    methods = {
        'VEP Plugin':  scored_df['nmd_vep_plugin'],
        'Lindeboom Rules': scored_df['nmd_rules'],
        'Consensus':   scored_df['nmd_consensus'],
    }
    x = range(len(cats))
    w = 0.25
    fig, ax = plt.subplots(figsize=(10, 4))
    for i, (name, series) in enumerate(methods.items()):
        counts = [int((series == c).sum()) for c in cats]
        offset = (i - 1) * w
        bars = ax.bar([xi + offset for xi in x], counts, w, label=name, alpha=0.85)
        for bar, n in zip(bars, counts):
            if n > 0:
                ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.05,
                        str(n), ha='center', fontsize=8)
    ax.set_xticks(list(x)); ax.set_xticklabels(cats, fontsize=9)
    ax.set_ylabel("Candidates"); ax.set_title("NMD classification per method", fontsize=12)
    ax.legend(fontsize=9); ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    return _b64fig(fig)


def plot_confidence(scored_df: pd.DataFrame) -> str:
    """Bar: confidence score distribution (0-3)."""
    scores = [3, 2, 1, 0]
    labels = ['3 — both agree', '2 — single method', '1 — disagree', '0 — no data']
    colors = ['#4caf50', '#ff9800', '#f44336', '#9e9e9e']
    counts = [int((scored_df['nmd_confidence_score'] == s).sum()) for s in scores]
    fig, ax = plt.subplots(figsize=(7, 3))
    bars = ax.bar(labels, counts, color=colors, width=0.5)
    for bar, n in zip(bars, counts):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.05,
                str(n), ha='center', fontsize=11, fontweight='bold')
    ax.set_ylabel("Candidates"); ax.set_title("Ensemble confidence score distribution", fontsize=12)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    return _b64fig(fig)

# ═════════════════════════════════════════════════════════════════════════════
# HTML REPORT
# ═════════════════════════════════════════════════════════════════════════════

def generate_report(scored_df: pd.DataFrame, out_dir: Path):
    # Generate plots
    img_tiers  = plot_tiers(scored_df)
    img_ic50   = plot_ic50(scored_df)
    img_method = plot_nmd_breakdown(scored_df)
    img_conf   = plot_confidence(scored_df)

    n_tot  = len(scored_df)
    n_s    = (scored_df['nmd_consensus'] == 'SENSITIVE').sum()
    n_i    = (scored_df['nmd_consensus'] == 'INSENSITIVE').sum()
    n_u    = scored_df['nmd_consensus'].isin(['UNCERTAIN','UNKNOWN']).sum()
    n_t1   = (scored_df['priority_tier'] == 'TIER1').sum()
    n_t2   = (scored_df['priority_tier'] == 'TIER2').sum()
    n_dis  = (scored_df['nmd_confidence'] == 'methods_disagree').sum()
    hla_df = hla_allele_breakdown(scored_df)

    t1 = scored_df[scored_df['priority_tier'] == 'TIER1']
    t2 = scored_df[scored_df['priority_tier'] == 'TIER2']
    t3 = scored_df[scored_df['priority_tier'] == 'TIER3_control']

    COLS = ['Gene Name','HGVSp','Variant Type','hla_allele','MT Epitope Seq',
            'best_ic50','best_ic50_method','nmd_consensus','nmd_confidence',
            'nmd_confidence_score','nmd_vep_plugin','nmd_rules','nmd_rule_explanation']

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
           ".card{background:#f8f8f8;border-radius:8px;padding:14px 18px;display:inline-block;min-width:130px;margin:4px}"
           ".cv{font-size:22px;font-weight:500} .cl{font-size:11px;color:#777}"
           ".note{background:#fffbe6;border-left:3px solid #f0c040;padding:10px 14px;margin:12px 0;font-size:13px;color:#555}")

    def card(label, val, color=""):
        return (f'<div class="card"><div class="cl">{label}</div>'
                f'<div class="cv"{"" if not color else f" style=color:{color}"}>{val}</div></div>')

    def tbl(df, cols, empty="No candidates."):
        if df.empty:
            return f"<p style='color:#888;font-style:italic;'>{empty}</p>"
        th = ''.join(f'<th>{c}</th>' for c in cols)
        tr = ''.join('<tr>' + ''.join(f'<td>{r.get(c,"")}</td>' for c in cols) + '</tr>'
                     for _, r in df.iterrows())
        return f"<table><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table>"

    def hla_tbl():
        if hla_df.empty:
            return "<p style='color:#888;font-style:italic;'>No Tier 1 or Tier 2 candidates.</p>"
        th = ''.join(f'<th>{c}</th>' for c in ['HLA Allele','Candidates','Tier 1','Best IC50 (nM)','Median IC50 (nM)'])
        tr = ''.join(f'<tr><td>{r.hla_allele}</td><td>{r.n_candidates}</td><td>{r.n_tier1}</td>'
                     f'<td>{r.best_ic50:.1f}</td><td>{r.median_ic50:.1f}</td></tr>'
                     for _, r in hla_df.iterrows())
        return f"<table><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table>"

    def conf_cards():
        return ('<div style="display:flex;gap:12px;flex-wrap:wrap;margin-bottom:16px;">' +
                ''.join(f'<div class="card"><div class="cl">Score {s}</div>'
                        f'<div class="cv">{(scored_df.nmd_confidence_score==s).sum()}</div></div>'
                        for s in [3,2,1,0]) + '</div>')

    # Note: CSS braces are doubled {{ }} in f-string, image vars use .format() after
    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
<title>GBM NMD Scoring Report</title><style>{CSS}</style></head>
<body><div class="page">

<h1>GBM Pipeline — Stage 7: NMD Sensitivity Scoring</h1>
<p style="color:#888;font-size:13px;margin-top:4px;">
  GBM NMD-Neoantigen Pipeline &middot; github.com/paleslui/gbm-nmd-pipeline
</p>

<h2>What is this report?</h2>
<p class="sub">
  Each pVACseq neoantigen candidate is classified by NMD sensitivity using an ensemble of two
  methods: the VEP NMD plugin and Lindeboom et al. 2019 rule-based scoring.
  The core hypothesis is that TMZ-induced frameshift mutations create PTCs silenced by NMD —
  NMD inhibition could expose these neoantigens to immune recognition.
</p>

<h2>Summary</h2>
<div style="margin:16px 0;"><img src="{img_tiers}" style="max-width:85%;display:block;margin:0 auto;border-radius:4px;box-shadow:0 1px 4px rgba(0,0,0,.12);">
<p style="font-size:11px;color:#888;margin-top:6px;">Fig 1. Candidate counts per priority tier.</p></div>
{card("Total candidates", n_tot)}
{card("NMD-sensitive", n_s, COL_S)}
{card("NMD-insensitive", n_i, COL_I)}
{card("Uncertain / N/A", n_u, COL_U)}
{card("Tier 1", n_t1)}
{card("Tier 2", n_t2)}
{card("Methods disagree", n_dis)}

<h2>IC50 Distribution by NMD Class</h2>
<div style="margin:16px 0;"><img src="{img_ic50}" style="max-width:85%;display:block;margin:0 auto;border-radius:4px;box-shadow:0 1px 4px rgba(0,0,0,.12);">
<p style="font-size:11px;color:#888;margin-top:6px;">Fig 2. IC50 per candidate colored by NMD consensus. Dashed lines = Tier 1 (50nM) and Tier 2/3 (500nM) thresholds.</p></div>

<h2>Per-HLA Allele Breakdown</h2>
<p class="sub">Tier 1 and Tier 2 candidates across HLA alleles.
Alleles with multiple high-confidence binders are the strongest therapeutic targets.</p>
{hla_tbl()}

<h2>NMD Classification per Method</h2>
<div style="margin:16px 0;"><img src="{img_method}" style="max-width:85%;display:block;margin:0 auto;border-radius:4px;box-shadow:0 1px 4px rgba(0,0,0,.12);">
<p style="font-size:11px;color:#888;margin-top:6px;">Fig 3. NMD classification per scoring method — VEP plugin, Lindeboom rules, and ensemble consensus.</p></div>

<h2>Confidence Score Distribution</h2>
<div style="margin:16px 0;"><img src="{img_conf}" style="max-width:85%;display:block;margin:0 auto;border-radius:4px;box-shadow:0 1px 4px rgba(0,0,0,.12);">
<p style="font-size:11px;color:#888;margin-top:6px;">Fig 4. Ensemble confidence score distribution across all candidates.</p></div>
<p class="sub">Score 3 = both methods agree; 2 = single method; 1 = disagree; 0 = no data.</p>
{conf_cards()}

<h2>NMD Scoring Methods</h2>
<h3>Method 1 — VEP NMD Plugin</h3>
<p class="sub">
  The VEP NMD plugin (Ensembl v105+) annotates truncating variants with NMD_escaping_variant
  when the PTC is predicted to escape NMD. An empty NMD field on a truncating variant means
  NMD is triggered (SENSITIVE). Requires Dragen VEP-annotated VCF with --plugin NMD.
</p>

<h3>Method 2 — Lindeboom Rules (Nat Cell Biol 2019)</h3>
<p class="sub">Applied to frameshift and stop-gained variants. Rules in priority order:</p>
<ul style="font-size:13px;color:#444;margin:0 0 12px 20px;line-height:1.8">
  <li><strong>Rule 4 — Start-proximal PTC (&lt;150nt):</strong> Pioneer round completes before NMD surveillance → INSENSITIVE</li>
  <li><strong>Rule 1 — Last exon:</strong> No downstream EJC to trigger NMD → INSENSITIVE</li>
  <li><strong>Rule 3 — Long exon (&gt;407nt):</strong> EJC too far downstream → INSENSITIVE</li>
  <li><strong>Rule 2 — 55nt boundary:</strong> PTC &gt;55nt upstream of last EJC → SENSITIVE (canonical NMD)</li>
</ul>
<div class="note">⚠ NMD rules only apply to truncating mutations.
Missense variants (NOT_APPLICABLE) are not subject to NMD but remain immunogenic neoantigen candidates.</div>

<h2>Priority Tiers</h2>
<ul style="font-size:13px;color:#444;margin:0 0 12px 20px;line-height:1.8">
  <li><strong>Tier 1 — NMD-sensitive + IC50 &lt;50 nM:</strong> Primary therapeutic targets — NMD inhibition could expose these neoantigens.</li>
  <li><strong>Tier 2 — NMD-sensitive + IC50 50–500 nM:</strong> Moderate binders, potentially relevant after NMD inhibition.</li>
  <li><strong>Tier 3 — NMD-insensitive + IC50 &lt;500 nM:</strong> Already translated — controls for immune response without NMD inhibition.</li>
  <li><strong>Unclassified:</strong> Missense variants or insufficient transcript information.</li>
</ul>

<h2>Tier 1 — NMD-sensitive + IC50 &lt;50 nM</h2>
{tbl(t1, COLS, "No Tier 1 candidates in this sample.")}

<h2>Tier 2 — NMD-sensitive + IC50 50–500 nM</h2>
{tbl(t2, COLS, "No Tier 2 candidates in this sample.")}

<h2>Tier 3 — NMD-insensitive controls</h2>
{tbl(t3, COLS, "No Tier 3 candidates in this sample.")}

<h2>All Candidates</h2>
<p class="sub">Complete table including missense variants (NOT_APPLICABLE).</p>
{tbl(scored_df, COLS)}

</div></body></html>"""

    out = out_dir / "report_nmd.html"
    out.write_text(html, encoding='utf-8')
    print(f"[REPORT] Saved {out}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="GBM Stage 7 — NMD Sensitivity Scoring")
    parser.add_argument('--pvacseq_tsv', required=True)
    parser.add_argument('--vep_vcf', default=None,
                        help='VEP-annotated VCF (gzipped). Optional — rule-based only if absent.')
    parser.add_argument('--out_dir', default='./nmd_output')
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}\n  GBM Pipeline - Stage 7: NMD Sensitivity Scoring\n{'='*60}\n")
    print(f"[INFO] Loading pVACseq results: {args.pvacseq_tsv}")
    pvacseq_df = pd.read_csv(args.pvacseq_tsv, sep='\t')
    print(f"[INFO] {len(pvacseq_df)} candidates loaded")

    vep_variants = {}
    if args.vep_vcf:
        print(f"\n[INFO] Parsing VEP VCF: {args.vep_vcf}")
        vep_variants = parse_vep_vcf(Path(args.vep_vcf))
    else:
        print("\n[INFO] No VEP VCF provided — rule-based scoring only")

    print("\n[INFO] Scoring NMD sensitivity...")
    scored_df = score_candidates(pvacseq_df, vep_variants)
    scored_df.to_csv(out_dir / "nmd_scored_candidates.tsv", sep='\t', index=False)
    print(f"[TABLE] Saved {out_dir}/nmd_scored_candidates.tsv")

    hla_bd = hla_allele_breakdown(scored_df)
    if not hla_bd.empty:
        hla_bd.to_csv(out_dir / "nmd_hla_breakdown.tsv", sep='\t', index=False)
        print(f"[TABLE] Saved {out_dir}/nmd_hla_breakdown.tsv")

    print(f"\n{'─'*40}")
    for label, val in [
        ("Total candidates",   len(scored_df)),
        ("NMD-SENSITIVE",      (scored_df['nmd_consensus']=='SENSITIVE').sum()),
        ("NMD-INSENSITIVE",    (scored_df['nmd_consensus']=='INSENSITIVE').sum()),
        ("NOT APPLICABLE",     (scored_df['nmd_rules']=='NOT_APPLICABLE').sum()),
        ("UNCERTAIN/UNKNOWN",  scored_df['nmd_consensus'].isin(['UNCERTAIN','UNKNOWN']).sum()),
        ("Methods disagree",   (scored_df['nmd_confidence']=='methods_disagree').sum()),
        ("Tier 1",             (scored_df['priority_tier']=='TIER1').sum()),
        ("Tier 2",             (scored_df['priority_tier']=='TIER2').sum()),
        ("Tier 3 controls",    (scored_df['priority_tier']=='TIER3_control').sum()),
    ]:
        print(f"  {label:<22} {val}")
    print(f"{'─'*40}")

    print("\n[REPORT] Generating HTML report...")
    generate_report(scored_df, out_dir)
    print(f"\n{'='*60}\n  Stage 7 complete. Output: {out_dir.resolve()}\n{'='*60}\n")


if __name__ == '__main__':
    main()