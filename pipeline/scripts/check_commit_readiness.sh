#!/usr/bin/env bash
# =============================================================================
# check_commit_readiness.sh
# -----------------------------------------------------------------------------
# Diagnoses the working tree before a commit and tells you WHAT to commit.
#
# It separates two kinds of modified files:
#   1. REAL source changes  — logical code/content edits that belong in a commit
#   2. PATH-ONLY changes     — files that differ from HEAD only because setup.sh
#                              STEP 8 substituted the install path (and the
#                              generated conda.cacheDir line); safe to ignore /
#                              regenerate, and reset by reset_paths_for_commit.sh
#
# Classification method: for each of the 8 path-substituted files, normalise the
# working copy the same way reset_paths_for_commit.sh would (install path back to
# __BASE__, config.sh BASE= back to the placeholder, drop the setup-generated
# conda.cacheDir line) and diff against HEAD ignoring blank lines. If anything
# survives, the file has REAL changes on top of the path substitution.
#
# Usage:
#   bash pipeline/scripts/check_commit_readiness.sh
# =============================================================================

set -euo pipefail

PLACEHOLDER="/path/to/gbm_nmd_pipeline"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG="$REPO_ROOT/config.sh"

if ! git -C "$REPO_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "[ERROR] $REPO_ROOT is not a git work tree." >&2
    exit 1
fi

BASE="$(grep -m1 '^BASE=' "$CONFIG" | sed -e 's/^BASE=//' -e 's/^["'"'"']//' -e 's/["'"'"']$//')"

# The 8 files that carry per-install path substitutions
PATH_SUB_FILES=(
    "config.sh"
    "pipeline/nextflow-pvacseq/slurm.config"
    "pipeline/slurm/master_pipeline.sh"
    "pipeline/slurm/slurm_nmd.sh"
    "pipeline/slurm/slurm_pvacseq.sh"
    "pipeline/slurm/slurm_pvacseq_filtered.sh"
    "pipeline/slurm/slurm_python.sh"
    "pipeline/slurm/slurm_setup.sh"
)

is_path_sub_file() {
    local needle="$1" x
    for x in "${PATH_SUB_FILES[@]}"; do
        [[ "$x" == "$needle" ]] && return 0
    done
    return 1
}

# Produce the path-normalised working copy of a path-sub file on stdout
normalise() {
    local rel="$1"
    local out
    out="$(sed "s|${BASE}|__BASE__|g" "$REPO_ROOT/$rel")"
    if [[ "$rel" == "config.sh" ]]; then
        out="$(printf '%s\n' "$out" | sed "s|^BASE=.*|BASE=${PLACEHOLDER}|")"
    fi
    if [[ "$rel" == *"slurm.config" ]]; then
        out="$(printf '%s\n' "$out" | sed '/^conda\.cacheDir/d')"
    fi
    printf '%s\n' "$out"
}

REAL_FILES=()        # tracked, real logical changes
PATHONLY_FILES=()    # tracked, only path/setup-artifact changes
NEW_FILES=()         # untracked new files
PATHSUB_IN_DIFF=()   # path-sub files currently showing in the diff

# --- tracked modifications --------------------------------------------------
while IFS= read -r rel; do
    [[ -z "$rel" ]] && continue
    if is_path_sub_file "$rel"; then
        PATHSUB_IN_DIFF+=("$rel")
        committed="$(git -C "$REPO_ROOT" show "HEAD:$rel" 2>/dev/null || true)"
        if diff -bB <(printf '%s\n' "$committed") <(normalise "$rel") >/dev/null 2>&1; then
            PATHONLY_FILES+=("$rel")
        else
            REAL_FILES+=("$rel")
        fi
    else
        REAL_FILES+=("$rel")
    fi
done < <(git -C "$REPO_ROOT" diff --name-only)

# --- untracked new files ----------------------------------------------------
while IFS= read -r rel; do
    [[ -z "$rel" ]] && continue
    NEW_FILES+=("$rel")
done < <(git -C "$REPO_ROOT" ls-files --others --exclude-standard)

# --- report -----------------------------------------------------------------
echo "================================================================"
echo " Commit readiness — $REPO_ROOT"
echo " BASE in config.sh: $BASE"
echo "================================================================"
echo ""

echo "[1] REAL source changes (these belong in the commit):"
if [[ ${#REAL_FILES[@]} -eq 0 && ${#NEW_FILES[@]} -eq 0 ]]; then
    echo "    (none)"
else
    for f in "${REAL_FILES[@]}"; do echo "    M  $f"; done
    for f in "${NEW_FILES[@]}";  do echo "    A  $f  (new file)"; done
fi
echo ""

echo "[2] PATH-substitution only (safe to ignore — setup.sh STEP 8 artifacts):"
if [[ ${#PATHONLY_FILES[@]} -eq 0 ]]; then
    echo "    (none)"
else
    for f in "${PATHONLY_FILES[@]}"; do echo "    ~  $f"; done
fi
echo ""

echo "[3] Suggested commit command (source changes only):"
ADD_LIST=("${REAL_FILES[@]}" "${NEW_FILES[@]}")
if [[ ${#ADD_LIST[@]} -eq 0 ]]; then
    echo "    (nothing to add)"
else
    printf '    git add'
    for f in "${ADD_LIST[@]}"; do printf ' \\\n        %s' "$f"; done
    printf '\n'
fi
echo ""

echo "[4] Reminder:"
if [[ ${#PATHSUB_IN_DIFF[@]} -gt 0 && "$BASE" != "$PLACEHOLDER" ]]; then
    echo "    The following path-substituted file(s) still carry your install path"
    echo "    and appear in 'git diff':"
    for f in "${PATHSUB_IN_DIFF[@]}"; do echo "        $f"; done
    echo ""
    echo "    Run this FIRST to revert them to __BASE__ placeholders:"
    echo "        bash pipeline/scripts/reset_paths_for_commit.sh"
    echo ""
    echo "    Files in list [1] that are ALSO path-substituted (e.g. config.sh,"
    echo "    slurm_nmd.sh) carry real edits AND your install path — reset cleans"
    echo "    the path; your real edits remain staged for commit."
else
    echo "    Paths look clean (BASE == placeholder). No reset needed."
fi
echo ""
