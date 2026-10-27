# t02_candidate_checks.py -- track re-run, step 2: per-candidate annotation and consistency checks. No test.
# For every pVACseq candidate of run_20260701_085805 (the graded showcase run):
#  (1) fetch the VEP CSQ entry of the candidate's OWN transcript (pVACseq 'Transcript') and of the transcript the
#      original stage 3 used (canonical protein_coding, else first entry), with the original coordinate lookup
#      ((chrom,pos,ref,alt), else (chrom,pos+1,ref,alt)).
#  (2) check Ensembl 113 models against VEP: exon total, exon number of the variant (from cDNA_position) and the
#      5'UTR offset (cDNA_position - CDS_position).
#  (3) frameshift stop codon by three independent routes, as the CDS position of the first nt of the new stop:
#      A  HGVSp 'p.<aa><pos><aa>fsTer<N>'  -> stop aa = pos + N - 1
#      B  len(FrameshiftSequence) = L      -> stop aa = L + 1   (plugin sequence excludes the stop)
#      C  GRCh38 FASTA + Ensembl 113 exons, apply the VCF allele, read in frame to the first TAA/TAG/TGA
#     CDS position = 3*(aa-1)+1. All three must agree before a candidate is re-classified.
import gzip, re, glob, os, pickle
import numpy as np, pandas as pd

W = '/cfs/earth/scratch/paleslui/audit_scratch/track_rerun/'
R = '/cfs/earth/scratch/paleslui/gbm_nmd_pipeline/results/run_20260701_085805/2_pvacseq/'
FA = '/cfs/earth/scratch/paleslui/gbm_nmd_pipeline/data/reference/GRCh38.primary_assembly.fa'
M = pickle.load(open(W + 'tx_models_grch38.pkl', 'rb'))
TX = M['tx']

fai = {}
for l in open(FA + '.fai'):
    n, ln, off, lb, lbb = l.split('\t')[:5]
    fai[n] = (int(ln), int(off), int(lb), int(lbb))


def fetch(chrom, a, b):
    """1-based closed genomic interval from the indexed FASTA (FASTA names have no 'chr' prefix)."""
    ln, off, lb, lbb = fai[chrom if chrom in fai else chrom.replace('chr', '').replace('M', 'MT')]
    def o(p):
        return off + (p - 1) // lb * lbb + (p - 1) % lb
    with open(FA, 'rb') as fh:
        fh.seek(o(a))
        raw = fh.read(o(b) - o(a) + 1)
    return raw.replace(b'\n', b'').replace(b'\r', b'').decode().upper()


COMP = str.maketrans('ACGTN', 'TGCAN')
STOPS = {'TAA', 'TAG', 'TGA'}


def route_c(t, chrom, pos, ref, alt):
    """Stop aa index of the mutated CDS (1-based), reading from the start codon."""
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    ex = t['exons']
    pieces = []
    for a, b in ex:
        s = fetch(chrom, a, b)
        if a <= pos <= b + 1 and (ref == '' and a <= pos - 1 <= b or ref != '' and pos + len(ref) - 1 <= b):
            if ref and s[pos - a:pos - a + len(ref)] != ref:
                return None, 'refmismatch'
            s = s[:pos - a] + alt + s[pos - a + len(ref):]
        pieces.append(s if t['strand'] == '+' else s.translate(COMP)[::-1])
    seq = ''.join(pieces)
    cs = t['cds_start_tx']
    for i in range(cs - 1, len(seq) - 2, 3):
        if seq[i:i + 3] in STOPS:
            return (i - (cs - 1)) // 3 + 1, 'ok'
    return None, 'nonstop'


rows = []
for f in sorted(glob.glob(R + 'pvactools/*/MHC_Class_I/*.filtered.tsv')):
    sample = os.path.basename(f).replace('.filtered.tsv', '')
    try:
        d = pd.read_csv(f, sep='\t', low_memory=False)
    except pd.errors.EmptyDataError:
        continue
    if d.empty:
        continue
    vcf = R + f'ensemblvep/{sample.replace("_", "")}-ensemble-annotated_vep.vcf.gz'
    fields, csq = None, {}
    for line in gzip.open(vcf, 'rt'):
        if line.startswith('##INFO=<ID=CSQ'):
            fields = re.search(r'Format: ([^"]+)"', line).group(1).strip().split('|')
        if line.startswith('#'):
            continue
        p = line.strip().split('\t')
        c = next((x[4:] for x in p[7].split(';') if x.startswith('CSQ=')), '')
        csq[(p[0], int(p[1]), p[3], p[4])] = [dict(zip(fields, e.split('|'))) for e in c.split(',')
                                              if len(e.split('|')) == len(fields)]
    for i, r in d.iterrows():
        key = (r.Chromosome, int(r.Start), r.Reference, r.Variant)
        tl = csq.get(key)
        vkey = key
        if tl is None:
            vkey = (r.Chromosome, int(r.Start) + 1, r.Reference, r.Variant)
            tl = csq.get(vkey, [])
        canon = [t for t in tl if t.get('CANONICAL') == 'YES' and t.get('BIOTYPE') == 'protein_coding']
        used = canon[0] if canon else (tl[0] if tl else {})
        own = next((t for t in tl if t.get('Feature') == r.Transcript), {})
        rows.append(dict(sample=sample, row=i, chrom=r.Chromosome, start=int(r.Start), ref=r.Reference,
                         alt=r.Variant, vcf_pos=vkey[1] if tl else None, vtype=r['Variant Type'],
                         gene=r['Gene Name'], pv_tx=r.Transcript, orig_tx=used.get('Feature'),
                         own_found=bool(own), own_csq=own.get('Consequence'), own_exon=own.get('EXON'),
                         own_cdna=own.get('cDNA_position'), own_cds=own.get('CDS_position'),
                         own_hgvsp=own.get('HGVSp'), own_fsseq_len=len(own.get('FrameshiftSequence', '') or ''),
                         own_nmd=own.get('NMD'), orig_csq=used.get('Consequence'), orig_nmd=used.get('NMD')))
C = pd.DataFrame(rows)
print('candidates', len(C), ' own transcript found in VEP CSQ:', C.own_found.sum())


def first_int(s):
    m = re.match(r'(\d+)', str(s))
    return int(m.group(1)) if m else None


# ---- (2) model vs VEP consistency on the candidate's own transcript
chk = []
for r in C.itertuples():
    t = TX.get(r.pv_tx)
    if t is None or not r.own_found:
        chk.append((False, None, None, None)); continue
    ex_cur, ex_tot = (str(r.own_exon).split('/') + [None])[:2] if r.own_exon else (None, None)
    tot_ok = ex_tot is not None and int(ex_tot) == t['n_exons']
    cdna, cds = first_int(r.own_cdna), first_int(r.own_cds)
    off_ok = (cdna is not None and cds is not None and t['cds_start_tx'] is not None
              and cdna - cds == t['cds_start_tx'] - 1)
    cum = np.cumsum(t['exon_lens'])
    ex_of = int(np.searchsorted(cum, cdna) + 1) if cdna else None
    ex_ok = ex_cur is not None and ex_of == int(ex_cur)
    chk.append((True, tot_ok, off_ok, ex_ok))
C[['in_model', 'chk_exon_total', 'chk_utr_offset', 'chk_exon_number']] = pd.DataFrame(chk, index=C.index)
m = C[C.in_model]
print('in Ensembl 113 model:', len(m), '| exon total = VEP:', m.chk_exon_total.mean().round(4),
      '| 5UTR offset = VEP:', m.chk_utr_offset.mean().round(4), '| exon number = VEP:', m.chk_exon_number.mean().round(4))
print('transcript versions: complete CDS models', np.mean([TX[t]['complete_cds'] for t in m.pv_tx]).round(4))

# ---- (3) stop position, three routes, frameshift candidates whose own consequence is a frameshift
fs = C[C.own_csq.fillna('').str.contains('frameshift_variant')].copy()
A, B, Cc, st = [], [], [], []
for r in fs.itertuples():
    mm = re.search(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}fsTer(\d+)', str(r.own_hgvsp))
    A.append(int(mm.group(1)) + int(mm.group(2)) - 1 if mm else None)
    B.append(r.own_fsseq_len + 1 if r.own_fsseq_len else None)
    c, s = route_c(TX[r.pv_tx], r.chrom, int(r.vcf_pos), r.ref, r.alt) if r.pv_tx in TX and r.vcf_pos else (None, 'nomodel')
    Cc.append(c); st.append(s)
fs['stop_aa_A'], fs['stop_aa_B'], fs['stop_aa_C'], fs['route_c_status'] = A, B, Cc, st
fs['routes_agree'] = (fs.stop_aa_A == fs.stop_aa_B) & (fs.stop_aa_B == fs.stop_aa_C)
print(f'frameshift candidates on own transcript: {len(fs)} | A available {fs.stop_aa_A.notna().sum()} | '
      f'B available {fs.stop_aa_B.notna().sum()} | C ok {(fs.route_c_status == "ok").sum()} | all three agree {fs.routes_agree.sum()}')
print(fs[~fs.routes_agree][['sample', 'gene', 'pv_tx', 'own_hgvsp', 'stop_aa_A', 'stop_aa_B', 'stop_aa_C', 'route_c_status']].to_string(index=False))
C.to_csv(W + 'candidates_annot.tsv', sep='\t', index=False)
fs.to_csv(W + 'fs_stop_routes.tsv', sep='\t', index=False)
