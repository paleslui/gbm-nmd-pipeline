# Revision, October 2026: NMD labels at the PTC position

This folder records a correction to the NMD classification of the track cohort (graded run
`run_20260701_085805`, the run behind `results/showcase/`). The pipeline code and the showcase
results are left exactly as graded. The corrected labels and statistics are added here, next to them.

## What changed

**Old state.** `pipeline/nmd_scoring.py` placed the premature termination codon (PTC) at the
**position of the mutation**. For a nonsense (stop-gained) variant that is correct, because the
variant itself is the stop. For a **frameshift** it is not: the real PTC is the first stop codon
that the shifted reading frame meets downstream, often tens to hundreds of nucleotides later, and
sometimes in a different exon. Since the Lindeboom rules depend on where the PTC sits relative to
the exon junctions, frameshifts could get the wrong label.

**New state.** For frameshifts, the PTC is now placed at the **real stop codon**. Three independent
routes agree for all 99 frameshift candidates:

1. the VEP HGVSp `fsTer` offset
2. the length of the pVACseq frameshift sequence
3. a direct read of the mutated GRCh38 transcript (Ensembl 113) up to the first stop

Nonsense variants keep the variant position.

Two smaller fixes were made in the same pass:

- The rules are evaluated on the candidate's **own transcript** (the pVACseq `Transcript`), not
  on the canonical one.
- Rule 3 uses the length of the exon that actually holds the PTC, and rule 2 applies the 55-nt
  distance test explicitly.

The rule thresholds, the rule order, the tier logic and `nmd_cohort_summary.py` are unchanged.

## Labellings compared (cumulative)

| Label | Content |
|---|---|
| L0 | original (reproduced exactly from the graded run) |
| L1 | + own transcript |
| L2 | + PTC at the real stop for frameshifts |
| L3 | + rule 2/3 fixes, ensemble with the VEP plugin (secondary) |
| **L4** | **repaired Lindeboom rules alone, without the VEP plugin (primary)** |

The VEP plugin is dropped in L4 because it is Ensembl code that still works from the variant
position, so it cannot be repaired from our side.

## Effect

Labels of the 99 frameshift candidates, L0 to L4: 19 of 62 SENSITIVE became INSENSITIVE, 1
INSENSITIVE became SENSITIVE, 4 UNCERTAIN became INSENSITIVE, and 7 previously unlabelled ones got
a label. In total, 63 SENSITIVE / 26 INSENSITIVE became 48 SENSITIVE / 51 INSENSITIVE. Missense
and in-frame candidates are unaffected (one missense label returns to unclassified).

Paired primary vs recurrent tests (Wilcoxon signed-rank, BH over the 11 stage-3 metrics, exactly as
in `nmd_cohort_summary.py`):

| Metric | p original (L0) | p repaired (L4) | BH-adjusted, L4 |
|---|---|---|---|
| nmd_fraction | 0.063 | 0.168 | 0.231 |
| nmd_sensitive | 0.102 | 0.248 | 0.273 |
| fraction_escape_all | 0.031 | 0.168 | 0.231 |
| nmd_sensitive_neo_strong | 0.039 | 0.150 | 0.231 |
| fraction_escape_strong | 0.125 | 0.375 | 0.375 |
| Module 4 (neoepitope burden) | unchanged | unchanged | unchanged |

**Conclusion.** None of the label-dependent results was significant after BH correction in the
original run either. The two that were nominally significant (`fraction_escape_all`,
`nmd_sensitive_neo_strong`) are no longer nominally significant after the repair. The
neoepitope-burden results (Module 4: total peptides and strong binders, BH p = 0.024) do not depend
on the NMD label and are unchanged.

All labellings side by side are in `paired_stats_all_labellings.tsv`.

## Files

| Path | Content |
|---|---|
| `scripts/t01_tx_models_grch38.py` | builds Ensembl 113 GRCh38 transcript models |
| `scripts/t02_candidate_checks.py` | per-candidate annotation, three-route stop-codon check |
| `scripts/t03_relabel.py` | applies labellings L1 to L4 and reruns the unchanged cohort summary |
| `scripts/t04_compare.py` | side-by-side statistics of all labellings |
| `candidates_annot.tsv` | per-candidate annotation from t02 |
| `fs_stop_routes.tsv` | real stop position of each frameshift by the three routes |
| `labels_by_step.tsv` | label of every candidate under L0 to L4 |
| `paired_stats_all_labellings.tsv` | stage-3 paired statistics under every labelling |
| `results_L4/` | full stage-3 output (cohort and per sample) with the repaired labels, same layout as `results/showcase/3_nmd_analysis/` |

The scripts were run on the ZHAW cluster and contain absolute paths to that environment. They read
the graded run's stage-2 VEP/pVACseq output. Stages 1 and 2 were not rerun, since nothing in them
depends on the NMD label.

The same correction was applied in the follow-up GLASS analysis (master's thesis repository).
