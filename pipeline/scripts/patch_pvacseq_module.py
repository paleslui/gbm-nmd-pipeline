#!/usr/bin/env python3
"""
Patch the nf-core pVACseq Nextflow module to copy the IEDB binding-prediction
data directory to node-local scratch before invoking 'pvacseq run'.

WHY: NetMHCpan, NetMHCpanEL, and PickPocket each open many small data files
under iedb/mhc_i/ during prediction. On shared parallel filesystems (BeeGFS,
Lustre) under high concurrent I/O, transient open() failures can occur,
causing some predictions to silently fail and pVACseq's parser to crash with
'dict contains fields not in fieldnames' errors. Copying IEDB to node-local
scratch (e.g. /data/scratch/$SLURM_JOB_ID/) eliminates this race.

Same algorithms, same logic — purely a cluster-portability adjustment.
The patch is idempotent: running it twice is safe (it detects an already-
patched module and exits cleanly).

Usage:
  python3 patch_pvacseq_module.py [path_to_main.nf]

Default path (when no arg): pipeline/nextflow-pvacseq/modules/local/pvacseq/main.nf
relative to current working directory.
"""
import sys
import os

DEFAULT_PATH = 'pipeline/nextflow-pvacseq/modules/local/pvacseq/main.nf'
SENTINEL = '# Localize IEDB to node-local scratch'

path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_PATH
if not os.path.isfile(path):
    print(f'ERROR: file not found: {path}')
    sys.exit(1)

with open(path) as f:
    content = f.read()

# Idempotency check
if SENTINEL in content:
    print(f'[OK] {path} already patched (sentinel found)')
    sys.exit(0)

lines = content.splitlines(keepends=True)

# Find pvacseq run line
pvac_idx = None
for i, line in enumerate(lines):
    if line.startswith('    pvacseq run '):
        pvac_idx = i
        break
if pvac_idx is None:
    print('ERROR: could not find pvacseq run line in module')
    sys.exit(1)

# Find script body start (preceding line with bare triple-quote)
script_start = None
for i in range(pvac_idx - 1, -1, -1):
    if lines[i].strip() == '"""':
        script_start = i
        break
if script_start is None:
    print('ERROR: could not find script body start')
    sys.exit(1)

# Find blank line after pvacseq run block
pvac_end = None
for i in range(pvac_idx, len(lines)):
    if lines[i].strip() == '':
        pvac_end = i
        break
if pvac_end is None:
    print('ERROR: could not find pvacseq run block end')
    sys.exit(1)

# Replacement block (Nextflow Groovy escaping: \$ for bash vars; ${...} for Groovy interpolation)
new_text = '''    """
    # ----------------------------------------------------------------
    # Localize IEDB to node-local scratch to avoid shared-filesystem
    # file-open contention during NetMHCpan/NetMHCpanEL parallel
    # data file reads. Cluster-portability fix; same algorithms,
    # same logic. See pipeline/scripts/patch_pvacseq_module.py.
    # ----------------------------------------------------------------
    LOCAL_SCRATCH=\\${LSFM_CLUSTER_LOCAL_SCRATCH_ROOT_PATH:-/data/scratch}/\\${SLURM_JOB_ID:-\\$\\$}
    LOCAL_IEDB=\\$LOCAL_SCRATCH/iedb
    mkdir -p \\$LOCAL_SCRATCH
    if [ ! -d "\\$LOCAL_IEDB" ]; then
        echo "[\\$(date)] Copying IEDB to node-local scratch: \\$LOCAL_IEDB"
        cp -r ${iedb}/. \\$LOCAL_IEDB/
        echo "[\\$(date)] IEDB copy complete"
    fi

    pvacseq run \\\\
        ${vcf} \\\\
        ${sample_name} \\\\
        ${hla} \\\\
        ${algorithms} \\\\
        ${prefix} \\\\
        --iedb-install-directory \\$LOCAL_IEDB ${blastp_opt} ${genes_opt} ${phased_opt} ${args} -t ${task.cpus}

'''
# Note above: in this Python heredoc, \\$ becomes \$ in the file, which Nextflow's
# Groovy parser sees as \$ (literal backslash + dollar in Groovy → just dollar in bash).
# That keeps the var as a bash var, not Groovy interpolation.

new_lines = [line + '\n' for line in new_text.split('\n')[:-1]]
result = lines[:script_start] + new_lines + lines[pvac_end:]

with open(path, 'w') as f:
    f.writelines(result)

print(f'[OK] patched: {path}')
