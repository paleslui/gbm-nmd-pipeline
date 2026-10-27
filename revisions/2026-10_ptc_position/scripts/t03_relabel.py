# t03_relabel.py -- track re-run, step 3: repaired NMD labels for the graded track run (run_20260701_085805).
# =====================================================================================================
# PRE-SPECIFICATION (written before any repaired label or statistic was computed)
# Starting point: the reproduced original stage 3 (track_rerun/repro, identical to the graded run apart from a
#   dash character in one text column). Stages 1-2 are not re-run: nothing in them depends on the NMD label,
#   and every input the repair needs (candidate's own transcript, HGVSp, FrameshiftSequence, VEP NMD field) is
#   already in the stage-2 VEP files. The original nmd_cohort_summary.py is run UNCHANGED on every labelling.
# Labellings (cumulative; thresholds and rule order kept exactly as in the track: rule 4 (<=150 nt) > rule 1
#   (last exon) > rule 3 (>407 nt) > rule 2 (>55 nt upstream of the last junction); ensemble and tier logic
#   copied verbatim from nmd_scoring.py):
#   L0  original (reproduced)
#   L1  transcript fix: rules and VEP NMD plugin read on the candidate's OWN transcript (pVACseq 'Transcript')
#   L2  L1 + PTC at the real stop codon for frameshifts (all three stop routes of t02 agree for 99/99)
#   L3  L2 + rule 3 on the real length of the exon holding the PTC + rule 2 with the 55-nt distance test
#       -> repaired ensemble (VEP plugin unchanged: Ensembl code, still variant-position based) [SECONDARY]
#   L4  L3 rules alone, without the VEP plugin -> repaired rules [PRIMARY]
# Tests: the track's own label-dependent paired Wilcoxon signed-rank tests (T vs M, BH over the track's 11
#   stage-3 metrics, as computed by nmd_cohort_summary.py), for L4 (primary) and L3 (secondary):
#   Module 3: nmd_fraction, nmd_sensitive, nmd_insensitive; Module 5: fraction_escape_strong,
#   fraction_escape_all, nmd_sensitive_neo_strong, nmd_escape_neo_strong. 14 tests, logged to tests.tsv.
#   No expected direction is claimed; the question is whether the track's results survive the repair.
#   L1 and L2 are descriptive only (label changes and metric values, no tests).
# Least sure: the 2 frameshift candidates on 5'-incomplete Ensembl models (rule 4 is not meaningful there);
#   they are flagged, not dropped, to keep the candidate set identical to the original.
# =====================================================================================================
import re, os, pickle, shutil, glob
import numpy as np, pandas as pd

W = '/cfs/earth/scratch/paleslui/audit_scratch/track_rerun/'
TX = pickle.load(open(W + 'tx_models_grch38.pkl', 'rb'))['tx']
C = pd.read_csv(W + 'candidates_annot.tsv', sep='\t')
FS = pd.read_csv(W + 'fs_stop_routes.tsv', sep='\t').set_index(['sample', 'row'])
TRUNCATING = ['stop_gained', 'frameshift', 'stop_lost', 'start_lost', 'splice_donor', 'splice_acceptor']
RULE_EXPLANATIONS = {
    'rule1_last_exon':             'PTC in last exon - no downstream EJC',
    'rule2_55nt_boundary':         'PTC >55nt upstream of last exon junction - NMD-sensitive',
    'rule2_within_55nt':           'PTC within 55nt of last exon junction - escapes',
    'rule3_long_exon_>407nt':      'PTC in long exon (>407nt) - reduced EJC density',
    'rule4_start_proximal_<150nt': 'PTC within 150nt of start codon - ribosome reinitiation',
}


def first_int(s):
    m = re.match(r'(\d+)', str(s))
    return int(m.group(1)) if m else None


# ---------- verbatim from nmd_scoring.py (score_nmd_vep_plugin body, ensemble_nmd, tier rule)
def plugin(csq, nmd):
    if not any(c in str(csq) for c in TRUNCATING):
        return 'UNKNOWN'
    nmd = '' if (nmd is None or (isinstance(nmd, float) and np.isnan(nmd))) else str(nmd).strip()
    if not nmd:
        return 'SENSITIVE'
    return 'INSENSITIVE' if 'escaping' in nmd.lower() else 'UNKNOWN'


def ensemble_nmd(vep, rule):
    v = vep if vep in ('SENSITIVE', 'INSENSITIVE') else None
    r = rule if rule in ('SENSITIVE', 'INSENSITIVE') else None
    if v is None and r is None: return 'UNKNOWN', 'no_data', 0
    if v is None:               return r, 'single_method_rules', 2
    if r is None:               return v, 'single_method_vep', 2
    if v == r:                  return v, 'high_confidence', 3
    return 'UNCERTAIN', 'methods_disagree', 1


def tier(cons, ic50):
    return ('TIER1' if cons == 'SENSITIVE' and ic50 < 50 else 'TIER2' if cons == 'SENSITIVE' and ic50 < 500 else
            'TIER3_control' if cons == 'INSENSITIVE' and ic50 < 500 else 'UNCLASSIFIED')


# ---------- rules on the candidate's own transcript
def rules(r, fix_ptc, fix_rules):
    """Track rule order. Variant-position inputs exactly as the track read them (EXON, CDS_position,
    cDNA_position span) unless fixed."""
    if not any(c in str(r.own_csq) for c in TRUNCATING):
        return 'NOT_APPLICABLE', 'not_truncating', ''
    ex_cur, ex_tot = [int(x) for x in str(r.own_exon).split('/')] if '/' in str(r.own_exon) else (None, None)
    cds_var = first_int(r.own_cds)
    t = TX[r.pv_tx]
    if not fix_ptc and not fix_rules:
        cd = str(r.own_cdna).split('-')
        span = abs(int(cd[-1]) - int(cd[0])) + 1 if cd[0].isdigit() and cd[-1].isdigit() else None
        if cds_var is not None and cds_var <= 150: return 'INSENSITIVE', 'rule4_start_proximal_<150nt', ''
        if ex_cur is not None and ex_tot is not None:
            if ex_cur == ex_tot: return 'INSENSITIVE', 'rule1_last_exon', ''
            if span is not None and span > 407: return 'INSENSITIVE', 'rule3_long_exon_>407nt', ''
            if ex_tot - ex_cur >= 1: return 'SENSITIVE', 'rule2_55nt_boundary', ''
        return 'UNCERTAIN', 'insufficient_info', ''
    # PTC position (CDS coordinate of first nt of the stop codon, mutated transcript)
    if 'frameshift' in str(r.own_csq) and fix_ptc:
        stop_aa = FS.loc[(r.sample, r.row), 'stop_aa_C']
        p_cds = 3 * (int(stop_aa) - 1) + 1
    else:
        p_cds = 3 * ((cds_var - 1) // 3) + 1
    # mutated exon ends: exons from the variant's exon onwards shift by the indel length
    ref, alt = r.ref, r.alt
    while ref and alt and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
    delta = len(alt) - len(ref) if 'frameshift' in str(r.own_csq) else 0
    cum = np.cumsum(t['exon_lens'])
    k_v = int(np.searchsorted(cum, first_int(r.own_cdna)))          # 0-based exon of the variant
    cum_m = np.array([cum[k] + (delta if k >= k_v else 0) for k in range(len(cum))])
    p_tx = t['cds_start_tx'] + p_cds - 1
    k_p = int(np.searchsorted(cum_m, p_tx))                          # 0-based exon of the PTC
    n = len(cum_m)
    exon_len = cum_m[k_p] - (cum_m[k_p - 1] if k_p > 0 else 0)
    note = f'PTC CDS {p_cds}nt; exon {k_p + 1}/{n}; exon {exon_len}nt'
    if p_cds <= 150: return 'INSENSITIVE', 'rule4_start_proximal_<150nt', note
    if k_p == n - 1: return 'INSENSITIVE', 'rule1_last_exon', note
    if fix_rules:
        if exon_len > 407: return 'INSENSITIVE', 'rule3_long_exon_>407nt', note
        dist = cum_m[n - 2] - p_tx
        return ('SENSITIVE', 'rule2_55nt_boundary', note + f'; {dist}nt to last EJ') if dist > 55 else \
               ('INSENSITIVE', 'rule2_within_55nt', note + f'; {dist}nt to last EJ')
    # L2: real PTC exon, but the track's buggy rule 3 (variant cDNA span) and rule 2 (no distance test)
    cd = str(r.own_cdna).split('-')
    span = abs(int(cd[-1]) - int(cd[0])) + 1
    if span > 407: return 'INSENSITIVE', 'rule3_long_exon_>407nt', note
    return 'SENSITIVE', 'rule2_55nt_boundary', note


LAB = {'L1': (False, False, True), 'L2': (True, False, True), 'L3': (True, True, True), 'L4': (True, True, False)}
summary = []
for lab, (fix_ptc, fix_rules, use_plugin) in LAB.items():
    out = W + f'{lab}/per_sample'
    if os.path.exists(W + lab):
        shutil.rmtree(W + lab)
    os.makedirs(out)
    for sdir in sorted(glob.glob(W + 'repro/per_sample/*/')):
        s = os.path.basename(sdir.rstrip('/'))
        os.makedirs(f'{out}/{s}')
        src = sdir + 'nmd_scored_candidates.tsv'
        try:
            d = pd.read_csv(src, sep='\t', dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError:
            shutil.copyfile(src, f'{out}/{s}/nmd_scored_candidates.tsv')   # 0-candidate sample, kept as is
            continue
        cs = C[C['sample'] == s].sort_values('row')
        assert len(cs) == len(d) and (cs.pv_tx.values == d['Transcript'].values).all() and \
               (cs.start.astype(str).values == d['Start'].values).all(), s
        for i, r in enumerate(cs.itertuples()):
            rcls, rapplied, notes = rules(r, fix_ptc, fix_rules)
            vep = plugin(r.own_csq, r.own_nmd) if use_plugin else None
            cons, conf, score = ensemble_nmd(vep, rcls)
            ic50 = float(d.at[i, 'best_ic50'])
            d.at[i, 'nmd_vep_plugin'] = vep if use_plugin else 'NOT_USED'
            d.at[i, 'nmd_rules'] = rcls
            d.at[i, 'nmd_rule_applied'] = rapplied
            d.at[i, 'nmd_rule_explanation'] = RULE_EXPLANATIONS.get(rapplied, rapplied or 'no transcript data')
            d.at[i, 'nmd_rule_notes'] = notes
            d.at[i, 'nmd_consensus'] = cons
            d.at[i, 'nmd_confidence'] = conf
            d.at[i, 'nmd_confidence_score'] = str(score)
            d.at[i, 'priority_tier'] = tier(cons, ic50)
            summary.append(dict(label=lab, sample=s, row=r.row, gene=r.gene, vtype=r.vtype, own_tx=r.pv_tx,
                                orig_tx=r.orig_tx, consensus=cons, rule=rapplied, rule_class=rcls, plugin=vep,
                                tier=tier(cons, ic50), incomplete_model=not TX[r.pv_tx]['complete_cds']))
        d.to_csv(f'{out}/{s}/nmd_scored_candidates.tsv', sep='\t', index=False)

S = pd.DataFrame(summary)
orig = []
for sdir in sorted(glob.glob(W + 'repro/per_sample/*/')):
    src = sdir + 'nmd_scored_candidates.tsv'
    try:
        d = pd.read_csv(src, sep='\t', dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        continue
    s = os.path.basename(sdir.rstrip('/'))
    for i, r in d.iterrows():
        orig.append(dict(label='L0', sample=s, row=i, vtype=r['Variant Type'], consensus=r.nmd_consensus,
                         rule=r.nmd_rule_applied, tier=r.priority_tier))
S = pd.concat([pd.DataFrame(orig), S], ignore_index=True)
S.to_csv(W + 'labels_by_step.tsv', sep='\t', index=False)
fs = S[S.vtype == 'FS']
print('frameshift neoantigen candidates: consensus by labelling')
print(pd.crosstab(fs.label, fs.consensus).to_string())
print('\nall candidates: tier by labelling')
print(pd.crosstab(S.label, S.tier).to_string())
w = fs.pivot_table(index=['sample', 'row'], columns='label', values='consensus', aggfunc='first')
for a, b in (('L0', 'L1'), ('L1', 'L2'), ('L2', 'L3'), ('L0', 'L3'), ('L0', 'L4')):
    print(f'{a}->{b}: frameshift labels changed {(w[a] != w[b]).sum()} of {len(w)}')
print('\nL0 vs L4 frameshift labels:'); print(pd.crosstab(w.L0, w.L4).to_string())
print('\nL4 rule that decided each frameshift:'); print(fs[fs.label == 'L4'].rule.value_counts().to_string())
print('non-FS candidates whose label changed L0->L4:',
      (S[(S.vtype != 'FS')].pivot_table(index=['sample', 'row'], columns='label', values='consensus', aggfunc='first')
       .pipe(lambda x: (x.L0 != x.L4).sum())))
