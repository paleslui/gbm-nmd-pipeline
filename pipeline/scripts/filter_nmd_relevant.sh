#!/usr/bin/env bash
# filter_nmd_relevant.sh
#
# Filter VEP-annotated VCFs to keep only variants relevant for NMD analysis,
# dropping pure missense / synonymous / UTR / intron / intergenic variants.
#
# A variant is KEPT if ANY of its VEP CSQ entries (or SnpEff ANN entries)
# contains one of: frameshift, stop_gained, stop_lost, start_lost,
# splice_donor, splice_acceptor, inframe_insertion, inframe_deletion.
#
# A variant is also kept if it carries a splice_region_variant overlap
# alongside a missense (these can affect splicing → potential PTC).
#
# Usage:
#   filter_nmd_relevant.sh <input.vcf[.gz]> <output.vcf.gz>
#
# Requires: bcftools (for indexing); pure awk for filtering logic.
#
# Author: Luigi Palese, GBM NMD-Neoantigen Pipeline (v2)

set -euo pipefail

INPUT="${1:?Usage: $0 <input.vcf[.gz]> <output.vcf.gz>}"
OUTPUT="${2:?Usage: $0 <input.vcf[.gz]> <output.vcf.gz>}"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input file not found: $INPUT" >&2
    exit 1
fi

# Auto-detect compression
if [[ "$INPUT" == *.gz ]]; then
    READ_CMD="zcat"
else
    READ_CMD="cat"
fi

# Whether output should be gzipped (recommended)
if [[ "$OUTPUT" == *.gz ]]; then
    WRITE_CMD="bgzip -c"
else
    WRITE_CMD="cat"
fi

# Track stats
TMP_STATS=$(mktemp)
trap 'rm -f "$TMP_STATS"' EXIT

$READ_CMD "$INPUT" | awk -v stats="$TMP_STATS" '
BEGIN {
    # NMD-relevant consequence patterns. ANY match → keep variant.
    # Order matters only for documentation; awk tests them all.
    keep_patterns["frameshift_variant"]      = 1
    keep_patterns["stop_gained"]             = 1
    keep_patterns["stop_lost"]               = 1
    keep_patterns["start_lost"]              = 1
    keep_patterns["splice_donor_variant"]    = 1
    keep_patterns["splice_acceptor_variant"] = 1
    keep_patterns["inframe_insertion"]       = 1
    keep_patterns["inframe_deletion"]        = 1
    keep_patterns["protein_altering_variant"]= 1   # rare, indel-like

    total = 0; kept = 0
    counts["frameshift_variant"] = 0
    counts["stop_gained"] = 0
    counts["stop_lost"] = 0
    counts["start_lost"] = 0
    counts["splice_donor_variant"] = 0
    counts["splice_acceptor_variant"] = 0
    counts["inframe_insertion"] = 0
    counts["inframe_deletion"] = 0
    counts["protein_altering_variant"] = 0
    counts["splice_region_only"] = 0
}

# Pass through all header lines unchanged
/^#/ { print; next }

# Body: one variant per line
{
    total++
    line = $0
    keep = 0
    matched_term = ""

    # Test each NMD-relevant pattern. ANY match → keep.
    for (pat in keep_patterns) {
        if (line ~ pat) {
            keep = 1
            matched_term = pat
            counts[pat]++
            break
        }
    }

    # Edge case: splice_region_variant alongside missense — affects splicing,
    # could create PTC. Keep these too.
    if (!keep && line ~ /splice_region_variant/ && line ~ /missense_variant/) {
        keep = 1
        counts["splice_region_only"]++
    }

    if (keep) {
        print
        kept++
    }
}

END {
    printf "TOTAL_VARIANTS\t%d\n", total > stats
    printf "KEPT_VARIANTS\t%d\n",  kept  >> stats
    printf "DROPPED_VARIANTS\t%d\n", total - kept >> stats
    printf "RETAIN_RATE\t%.2f%%\n", 100.0 * kept / total >> stats
    printf "BREAKDOWN_BY_FIRST_MATCH\n" >> stats
    for (p in counts) {
        if (counts[p] > 0)
            printf "  %s\t%d\n", p, counts[p] >> stats
    }
}
' | $WRITE_CMD > "$OUTPUT"

# Index the output if gzipped
if [[ "$OUTPUT" == *.gz ]]; then
    bcftools index --tbi -f "$OUTPUT" 2>/dev/null || tabix -p vcf -f "$OUTPUT"
fi

# Print stats
echo ""
echo "=== Filter results: $(basename "$INPUT") ==="
cat "$TMP_STATS"
echo "Output: $OUTPUT"