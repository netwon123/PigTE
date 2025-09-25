#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# TE→VCF annotator
# -----------------------------
# Usage:
#   te_annotate.sh \
#     -t Best_TE_annotation_per_SV_ref-TE.tsv \
#     -v panel_del.vcf.gz \
#     -r /path/to/reference.fa.fai \
#     -o panel_del \
#     [-H header.txt]
#
# Outputs:
#   panel_del.annotated.vcf.gz
#   panel_del.annotated.vcf.gz.tbi
#   panel_del.TEonly.vcf.gz
#   panel_del.TEonly.vcf.gz.tbi
#
# Notes:
# - Requires: awk, bedtools, bgzip, tabix, bcftools
# - If -H not provided, a minimal header defining INFO/TE_ID and INFO/TE_CLASS is auto-generated.

usage() {
  echo "Usage: $0 -t TE_TSV -v VCF[.gz] -r REF.FAI -o OUTPREFIX [-H HEADER.txt]" 1>&2
  exit 1
}

TE_TSV=""
VCF_IN=""
FAI=""
OUTPREFIX=""
HEADER=""

while getopts ":t:v:r:o:H:" opt; do
  case "$opt" in
    t) TE_TSV="$OPTARG" ;;
    v) VCF_IN="$OPTARG" ;;
    r) FAI="$OPTARG" ;;
    o) OUTPREFIX="$OPTARG" ;;
    H) HEADER="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -z "$TE_TSV" || -z "$VCF_IN" || -z "$FAI" || -z "$OUTPREFIX" ]] && usage
[[ -f "$TE_TSV" ]] || { echo "ERROR: TE TSV not found: $TE_TSV" >&2; exit 2; }
[[ -f "$VCF_IN" ]] || { echo "ERROR: VCF not found: $VCF_IN" >&2; exit 2; }
[[ -f "$FAI" ]] || { echo "ERROR: .fai not found: $FAI" >&2; exit 2; }

# tool checks
for x in awk bedtools bgzip tabix bcftools; do
  command -v "$x" >/dev/null 2>&1 || { echo "ERROR: $x not found in PATH" >&2; exit 3; }
done

# temp dir
TMPDIR="$(mktemp -d -t teanno.XXXXXX)"
trap 'rm -rf "$TMPDIR"' EXIT

# 1) Build CHROM,POS,END,TE_ID,TE_CLASS table from the TE TSV
#    Expect col1 like "SV:chr-pos" (split on ':' and '-') and TE_ID in col4, TE_CLASS in col5.
awk -F'\t' 'NR==1{print "CHROM\tPOS\tEND\tTE_ID\tTE_CLASS"; next} {split($1,a,"[:-]"); print a[2]"\t"a[3]"\t"a[3]"\t"$4"\t"$5}' \
  "$TE_TSV" > "$TMPDIR/te.table"

# 2) Remove header and sort with faidx order
tail -n +2 "$TMPDIR/te.table" > "$TMPDIR/te.table.nohdr"
bedtools sort -faidx "$FAI" -i "$TMPDIR/te.table.nohdr" > "$TMPDIR/te.table.sorted"

# 3) Reduce to CHROM,POS,TE_ID,TE_CLASS for site-level annotation (POS used as start=end)
awk -v OFS='\t' '{print $1,$2,$4,$5}' "$TMPDIR/te.table.sorted" > "$TMPDIR/te.bed2"

# 4) Compress and index
bgzip -c "$TMPDIR/te.bed2" > "$TMPDIR/te.bed2.gz"
tabix -s1 -b2 -e2 "$TMPDIR/te.bed2.gz"

# 5) Prepare header (if not provided)
HDR="$TMPDIR/add.header.txt"
if [[ -n "$HEADER" ]]; then
  [[ -f "$HEADER" ]] || { echo "ERROR: header.txt not found: $HEADER" >&2; exit 4; }
  cp "$HEADER" "$HDR"
else
  cat > "$HDR" <<'EOF'
##INFO=<ID=TE_ID,Number=1,Type=String,Description="Transposable element ID from TE annotation">
##INFO=<ID=TE_CLASS,Number=1,Type=String,Description="Transposable element class from TE annotation">
EOF
fi

# 6) Annotate
bcftools annotate \
  -a "$TMPDIR/te.bed2.gz" \
  -c CHROM,POS,INFO/TE_ID,INFO/TE_CLASS \
  -h "$HDR" \
  "$VCF_IN" \
  -O z -o "${OUTPREFIX}.annotated.vcf.gz"

tabix -p vcf "${OUTPREFIX}.annotated.vcf.gz"

# 7) Extract TE-only sites
bcftools view -O z -o "${OUTPREFIX}.TEonly.vcf.gz" \
  -i 'INFO/TE_ID!=""' "${OUTPREFIX}.annotated.vcf.gz"

tabix -p vcf "${OUTPREFIX}.TEonly.vcf.gz"

echo "DONE."
echo "Annotated VCF : ${OUTPREFIX}.annotated.vcf.gz"
echo "TE-only VCF   : ${OUTPREFIX}.TEonly.vcf.gz"
