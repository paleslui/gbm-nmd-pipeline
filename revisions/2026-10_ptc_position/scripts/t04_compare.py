# t04_compare.py -- track re-run, step 4: side-by-side of the track's stage-3 paired statistics under each
# labelling, and logging of the 14 tests pre-specified in t03_relabel.py (L4 primary, L3 secondary; the 7
# label-dependent metrics each). p and BH-adjusted p are exactly as computed by the unchanged
# nmd_cohort_summary.py. Module 4 metrics do not depend on the label and must be identical in every labelling.
import sys
import pandas as pd
sys.path.insert(0, '/cfs/earth/scratch/paleslui/audit_scratch')
from testlog import log

W = '/cfs/earth/scratch/paleslui/audit_scratch/track_rerun/'
ORIG = '/cfs/earth/scratch/paleslui/gbm_nmd_pipeline/results/run_20260701_085805/3_nmd_analysis/cohort/'
LABEL_DEP = ['nmd_fraction', 'nmd_sensitive', 'nmd_insensitive', 'fraction_escape_strong', 'fraction_escape_all',
             'nmd_sensitive_neo_strong', 'nmd_escape_neo_strong']
tabs = {'L0 original (graded run)': pd.read_csv(ORIG + 'paired_stats_stage3.tsv', sep='\t')}
for lab, name in (('repro', 'L0 reproduced'), ('L1', 'L1 own transcript'), ('L2', 'L2 + real PTC'),
                  ('L3', 'L3 + rule fixes, with VEP plugin'), ('L4', 'L4 repaired rules only')):
    tabs[name] = pd.read_csv(W + f'{lab}/cohort/paired_stats_stage3.tsv', sep='\t')

rows = []
for name, t in tabs.items():
    for r in t.itertuples():
        rows.append(dict(labelling=name, module=r.module, metric=r.metric, n_pairs=r.n_pairs,
                         median_T=r.median_primary, median_M=r.median_recurrent, mean_delta=r.mean_delta,
                         p=r.p_value, p_adj=r.p_value_adj_bh))
A = pd.DataFrame(rows)
A.to_csv(W + 'paired_stats_all_labellings.tsv', sep='\t', index=False)
m4 = A[A.module == 4].pivot_table(index='metric', columns='labelling', values='p')
assert (m4.nunique(axis=1) == 1).all(), 'Module 4 changed - should not depend on the label'
print('Module 4 (label-independent) identical across labellings: OK')
pd.set_option('display.width', 250)
for metric in LABEL_DEP:
    x = A[A.metric == metric].set_index('labelling')[['n_pairs', 'median_T', 'median_M', 'mean_delta', 'p', 'p_adj']]
    print(f'\n--- {metric}')
    print(x.to_string(float_format=lambda v: f'{v:.4g}'))

for lab, tag, role in (('L4 repaired rules only', 'TR-L4', 'primary'), ('L3 + rule fixes, with VEP plugin', 'TR-L3', 'secondary')):
    t = A[A.labelling == lab].set_index('metric')
    for metric in LABEL_DEP:
        r = t.loc[metric]
        log(f'{tag}-{metric}', 't04_compare.py', 'track', f'track re-run ({role}): {metric}, T vs M',
            'track 28-pair cohort, graded run 20260701', 'Wilcoxon signed-rank (paired), BH over 11 stage-3 metrics',
            f'{int(r.n_pairs)} pairs', float(r.mean_delta), f'mean delta M-T; BH p = {r.p_adj:.3g}', float(r.p),
            'none claimed (does the track result survive)')
print('\nlogged 14 tests')
