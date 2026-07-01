# Troubleshooting

Common failure modes when running the pipeline, and how to diagnose them.

These are the reproducibility failure modes hit during level-3/level-4
debugging. Each entry lists the symptom, the root cause, and the concrete fix.

## Setup phase

### `[STEP 10 FAILED]` — env integrity gate refuses to pass

Step 10 of `slurm_setup.sh` imports a probe package from each conda env as a
final integrity gate. A failure means one of the conda envs is corrupted or
missing packages (usually a partial or non-deterministic env build).

- **What to check:** which specific import failed — Step 10 prints the failing
  package and the env it belongs to.
- **Root cause pattern:** conda non-determinism during env build (a solve that
  succeeds once and drops a package the next time).
- **Fix:** rebuild the affected env from scratch:
  ```bash
  rm -rf conda_envs/nf_pvacseq nextflow_work/conda/env-4a82*
  sbatch pipeline/slurm/slurm_setup.sh
  ```

### `ModuleNotFoundError` mentions scipy, Bio, or pvactools during Stage 1 or Stage 2

The cluster miniconda hijacks `python` on compute nodes.

- **Symptom:** works when tested interactively but fails when run via SLURM.
- **Root cause:** `module load lsfm-init-miniconda/1.0.0` prepends the cluster
  python 3.9 to `PATH`, so a bare `python` on a compute node is not the env
  python that has scipy / Bio / pvactools installed.
- **Fix:** verify `slurm_python.sh` and `slurm_nmd.sh` invoke
  `$CONDA_PREFIX/bin/python` (not bare `python`).
- **Verify:** the `[INFO] Using python:` line in the SLURM log should show
  `conda_envs/nf_pvacseq/bin/python`.

## Stage 2 (pVACseq) phase

### Many `.command.err` files containing `No module named 'Bio'`

This was the big one. The pVACseq conda env (`env-4a82*`) spawns an MHCnuggets
child process with a bare `python`, which resolves to the wrong interpreter and
cannot import Biopython.

- **Root cause:** `prediction_class.py` at line 156 spawns the MHCnuggets child
  with bare `"python"` instead of the current interpreter.
- **The patch that fixes it:** rewrites that call to use `sys.executable`, and is
  applied by `setup.sh` Step 4b.
- **If the patch is missing after a rebuild:** check that `setup.sh` Step 4b
  actually ran and that its idempotency guard did not silently skip re-patching
  a freshly rebuilt env.
- **Diagnostic:**
  ```bash
  grep "sys.executable" nextflow_work/conda/env-4a82*/lib/python*/site-packages/pvactools/lib/prediction_class.py
  ```
  No match means the patch is missing — re-run `setup.sh` Step 4b.

### All pVACseq tasks fail after ~5 minutes

Usually the MHCnuggets binding-prediction step failing across the board.

- **Check the Nextflow trace:**
  `results/run_*/2_pvacseq/pipeline_info/execution_trace_*.txt`
- Grep for `FAILED` lines and note the task hash.
- Look at `.command.err` in the corresponding `nextflow_work/<hash>*/`
  directory for the underlying error (often the `No module named 'Bio'` issue
  above).

## Stage 3 phase

### `KeyError: 'priority_tier'` when running nmd_cohort_summary.py

- **Root cause:** a sample has zero candidates, but downstream code assumed the
  scored-candidates columns always exist.
- **Fix:** verified in Step 4 — the empty-sample guard in `nmd_scoring.py`
  `main()` returns early and writes only an empty `nmd_scored_candidates.tsv`,
  matching the baseline behavior.
- **Affected samples in our cohort (7 total):** 5_T, 5_M, 6_T, 17_T, 18_T,
  36_T, 52_M.

### Stage 3 crashes with `No filtered TSV files found in .../2_pvacseq/pvactools`

- Means Stage 2 produced no output even though it reported `COMPLETED`.
- Nextflow's `errorStrategy = 'retry then ignore'` can swallow all failures, so
  a stage that fails on every task still exits successfully.
- **Check the trace file** (above) for the `FAILED` task count — it should be
  zero. A non-zero count with an empty output dir means every task was ignored.

## Disk / quota

### setup.sh fails with `Disk quota exceeded` on your home partition

- **Root cause:** Nextflow writes its cache to `~/.nextflow` by default, which
  lives on the (small, quota-limited) home partition.
- **Fix:** `config.sh` exports `NXF_HOME=$BASE/.nextflow` to redirect the cache
  to scratch. Verify this line is present and that `BASE` points to scratch (not
  home).

## After committing / pushing

### Pipeline won't run after a fresh git pull

- **Symptom:** the working tree has `__BASE__` placeholders instead of real
  absolute paths.
- **Fix:** run `sbatch pipeline/slurm/slurm_setup.sh` — Step 8 re-applies the
  path substitutions. It is idempotent and takes ~seconds.
