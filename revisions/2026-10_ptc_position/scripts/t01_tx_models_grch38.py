# t01_tx_models_grch38.py -- track re-run, step 1: GRCh38 transcript models from the Ensembl 113 GTF
# (= GENCODE 47, the genebuild of the VEP 113 cache the track used). No statistical test.
# Same construction as audit a01: exons in transcript order, CDS span in transcript coordinates including the
# stop codon, 1-based closed coordinates. Incomplete models (cds_start_NF/cds_end_NF/mRNA_*_NF) are flagged.
import gzip, re, pickle
import numpy as np, pandas as pd

GTF = '/cfs/earth/scratch/paleslui/audit_scratch/downloads/ensembl113/Homo_sapiens.GRCh38.113.gtf.gz'
OUT = '/cfs/earth/scratch/paleslui/audit_scratch/track_rerun/tx_models_grch38.pkl'
attr_re = re.compile(r'(\S+) "([^"]*)"')
tx = {}
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        feat = p[2]
        if feat not in ('transcript', 'exon', 'CDS', 'start_codon', 'stop_codon'):
            continue
        tid = re.search(r'transcript_id "([^"]+)"', p[8]).group(1)
        if feat == 'transcript':
            d, tags = {}, []
            for k, v in attr_re.findall(p[8]):
                (tags.append(v) if k == 'tag' else d.__setitem__(k, v))
            tx[tid] = dict(tid=tid, version=d.get('transcript_version'), gene_name=d.get('gene_name', ''),
                           biotype=d.get('transcript_biotype'), chrom='chr' + p[0] if p[0] != 'MT' else 'chrM',
                           strand=p[6], tags=tags, exons=[], cds=[], start=[], stop=[])
        else:
            key = {'exon': 'exons', 'CDS': 'cds', 'start_codon': 'start', 'stop_codon': 'stop'}[feat]
            tx[tid][key].append((int(p[3]), int(p[4])))

rows = []
for tid, t in tx.items():
    s = t['strand']
    t['exons'] = sorted(t['exons'], key=lambda e: e[0], reverse=(s == '-'))
    ex = t['exons']
    t['n_exons'] = len(ex)
    t['exon_lens'] = [b - a + 1 for a, b in ex]
    bad = {'cds_start_NF', 'cds_end_NF', 'mRNA_start_NF', 'mRNA_end_NF'}
    t['complete_cds'] = bool(t['cds']) and bool(t['start']) and bool(t['stop']) and not (bad & set(t['tags']))
    t['cds_start_tx'] = t['cds_end_tx'] = None
    if t['cds']:
        g = [x for a, b in t['cds'] + t['stop'] for x in range(a, b + 1)]
        def g2t(gp):
            off = 0
            for a, b in ex:
                if a <= gp <= b:
                    return off + (gp - a + 1 if s == '+' else b - gp + 1)
                off += b - a + 1
        lo, hi = (min(g), max(g)) if s == '+' else (max(g), min(g))
        t['cds_start_tx'], t['cds_end_tx'] = g2t(lo), g2t(hi)
    rows.append(dict(tid=tid, version=t['version'], gene=t['gene_name'], biotype=t['biotype'], chrom=t['chrom'],
                     strand=s, n_exons=t['n_exons'], complete_cds=t['complete_cds'],
                     cds_start_tx=t['cds_start_tx'], cds_end_tx=t['cds_end_tx']))
df = pd.DataFrame(rows)
pickle.dump({'tx': tx, 'df': df}, open(OUT, 'wb'))
print(df.shape, 'protein_coding', (df.biotype == 'protein_coding').sum(), 'complete CDS', df.complete_cds.sum())
