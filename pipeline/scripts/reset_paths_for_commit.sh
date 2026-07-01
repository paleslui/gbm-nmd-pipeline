#!/usr/bin/env bash
# =============================================================================
# reset_paths_for_commit.sh
# -----------------------------------------------------------------------------
# Reverts the per-install path substitutions that setup.sh STEP 8 bakes into the
# working tree, so the repo can be committed with clean __BASE__ placeholders.
#
# setup.sh STEP 8 replaces the literal "__BASE__" placeholder (and config.sh's
# BASE= line) with the absolute install path. This script undoes that: it reads
# the current BASE from config.sh and rewrites that absolute path back to
# __BASE__ in the 8 affected files, then resets config.sh's BASE= line to the
# generic placeholder.
#
# It is IDEMPOTENT (safe to run repeatedly) and SELF-AWARE (does nothing if the
# tree is already reset).
#
# Usage:
#   bash pipeline/scripts/reset_paths_for_commit.sh            # apply
#   bash pipeline/scripts/reset_paths_for_commit.sh --dry-run  # preview only
#
# After running, commit ONLY the real source changes — see the final message.
# To restore your local paths afterwards, re-run setup.sh (idempotent) or its
# STEP 8 path-substitution block.
# =============================================================================

set -euo pipefail

PLACEHOLDER="/path/to/gbm_nmd_pipeline"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" || "${1:-}" == "-n" ]]; then
    DRY_RUN=1
fi

# Resolve repo root from this script's location (pipeline/scripts/ -> repo root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG="$REPO_ROOT/config.sh"

if [[ ! -f "$CONFIG" ]]; then
    echo "[ERROR] config.sh not found at $CONFIG" >&2
    exit 1
fi

# Read BASE= from config.sh, strip surrounding single/double quotes
BASE="$(grep -m1 '^BASE=' "$CONFIG" | sed -e 's/^BASE=//' -e 's/^["'"'"']//' -e 's/["'"'"']$//')"

if [[ -z "$BASE" ]]; then
    echo "[ERROR] Could not read BASE= from $CONFIG" >&2
    exit 1
fi

echo "[INFO] Repo root: $REPO_ROOT"
echo "[INFO] Current BASE in config.sh: $BASE"

# Self-aware: nothing to do if already at the generic placeholder
if [[ "$BASE" == "$PLACEHOLDER" ]]; then
    echo "[OK] BASE is already '$PLACEHOLDER' — already reset, no changes needed."
    exit 0
fi

# The 8 files that carry per-install path substitutions
FILES=(
    "config.sh"
    "pipeline/nextflow-pvacseq/slurm.config"
    "pipeline/slurm/master_pipeline.sh"
    "pipeline/slurm/slurm_nmd.sh"
    "pipeline/slurm/slurm_pvacseq.sh"
    "pipeline/slurm/slurm_pvacseq_filtered.sh"
    "pipeline/slurm/slurm_python.sh"
    "pipeline/slurm/slurm_setup.sh"
)

if [[ $DRY_RUN -eq 1 ]]; then
    echo "[DRY-RUN] No files will be modified. Showing planned changes:"
fi

modified=()
for rel in "${FILES[@]}"; do
    f="$REPO_ROOT/$rel"
    if [[ ! -f "$f" ]]; then
        echo "[WARN] Skipping missing file: $rel" >&2
        continue
    fi

    # Count how many lines would change (path occurrences + config.sh BASE line)
    hits=0
    if grep -q -- "$BASE" "$f"; then
        hits=$(grep -c -- "$BASE" "$f" || true)
    fi
    extra_base=0
    if [[ "$rel" == "config.sh" ]] && grep -q '^BASE=' "$f"; then
        extra_base=1
    fi

    if [[ $hits -eq 0 && $extra_base -eq 0 ]]; then
        continue
    fi

    if [[ $DRY_RUN -eq 1 ]]; then
        echo "  - $rel  (${hits} path occurrence(s)$([[ $extra_base -eq 1 ]] && echo ' + BASE= line'))"
        modified+=("$rel")
        continue
    fi

    # 1) Rewrite the absolute install path back to the __BASE__ placeholder
    sed -i "s|${BASE}|__BASE__|g" "$f"

    # 2) config.sh keeps a literal BASE= line (not __BASE__); reset it explicitly
    if [[ "$rel" == "config.sh" ]]; then
        sed -i "s|^BASE=.*|BASE=${PLACEHOLDER}|" "$f"
    fi

    # 3) slurm.config carries a setup.sh STEP 8-generated conda.cacheDir line
    #    (prepended to the file, with trailing blank lines) on top of the
    #    committed content. Strip it so the working copy matches HEAD byte-for-
    #    byte — otherwise the file keeps showing in 'git diff' after a reset.
    if [[ "$rel" == *"slurm.config" ]]; then
        sed -i '/^conda\.cacheDir/d' "$f"
        sed -i '/./,$!d' "$f"   # drop the leading blank line(s) it left behind
    fi

    modified+=("$rel")
done

echo ""
if [[ ${#modified[@]} -eq 0 ]]; then
    echo "[OK] No path-substituted files needed changes."
else
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "[DRY-RUN] ${#modified[@]} file(s) WOULD be modified:"
    else
        echo "[OK] Reset ${#modified[@]} file(s) to __BASE__ placeholders:"
    fi
    for m in "${modified[@]}"; do
        echo "    $m"
    done
fi

# Show a diff stat for the affected files (best-effort; needs a git repo)
if git -C "$REPO_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo ""
    echo "[INFO] git diff --stat for these files:"
    git -C "$REPO_ROOT" diff --stat -- "${FILES[@]}" || true
fi

echo ""
if [[ $DRY_RUN -eq 1 ]]; then
    echo "[DRY-RUN] Nothing was written. Re-run without --dry-run to apply."
    exit 0
fi

cat <<'EOF'

Paths reset. Review with 'git diff', stage the source-code changes
(gbm_analysis.py, nmd_cohort_summary.py, slurm_nmd.sh — these contain REAL code
changes, NOT just path changes), and commit. Then run setup.sh's Step 8 again to
restore your local paths (or re-run setup.sh entirely, which is idempotent).
EOF
