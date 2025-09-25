#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# TE→VCF annotator (with SVTYPE filter)
# -----------------------------
# Usage:
#   te_annotate.sh \
#     -t Best_TE_annotation_per_SV_ref-TE.tsv \
#     -v panel.vcf.gz \
#     -r /path/to/reference.fa.fai \
#     -o output_prefix \
#     -s INS|DEL \
#     [-H header.txt]
#
# Outputs (under OUTPREFIX.*):
#   OUTPREFIX.SVTYPE.annotated.vcf.gz
#   OUTPREFIX.SVTYPE.annotated.vcf.gz.tbi
#   OUTPREFIX.SVTYPE.TEonly.vcf.gz
#   OUTPREFIX.SVTYPE.TEonly.vcf.gz.tbi
#
# Requires: awk, bedtools, bgzip, tabix, bcftools

usage() {
  echo "Usage: $0 -t TE_TSV -v VCF[.gz] -r REF.FAI -o OUTPREFIX -s INS|DEL [-H HEADER.txt]" 1>&2
  exit 1
}

TE_TSV=""
VCF_IN=""
FAI=""
OUTPREFIX=""
HEADER=""
SVT=""

while getopts ":t:v:r:o:H:s:" opt; do
  case "$opt" in
    t) TE_TSV="$OPTARG" ;;
    v) VCF_IN="$OPTARG" ;;
    r) FAI="$OPTARG" ;;
    o) OUTPREFIX="$OPTARG" ;;
    H) HEADER="$OPTARG" ;;
    s) SVT="$OPTARG" ;;
    *) usage ;;

  esac
done

# shellcheck disable=SC2086
[[ -z "${TE_TSV}" || -z "${VCF_IN}" || -z "${FAI}" || -z "${OUTPREFIX}" || -z "${SVT}" ]] && usage
[[ "$SVT" == "INS" || "$SVT" == "DEL" ]] || { echo "ERROR: -s must be INS or DEL"; exit 2; }

[[ -f "$TE_TSV" ]] || { echo "ERROR: TE TSV not found: $TE_TSV" >&2; exit 2; }
[[ -f "$VCF_IN" ]] || { echo "ERROR: VCF not found: $VCF_IN" >&2; exit 2; }
[[ -f "$FAI" ]] || { echo "ERROR: .fai not found: $FAI" >&2; exit 2; }

for x in awk bedtools bgzip tabix bcftools; do
  command -v "$x" >/dev/null 2>&1 || { echo "ERROR: $x not found in PATH" >&2; exit 3; }
done

TMPDIR="$(mktemp -d -t teanno.XXXXXX)"
trap 'rm -rf "$TMPDIR"' EXIT

# 0) 先按 SVTYPE 过滤 VCF（兼容常见写法：INFO/SVTYPE=DEL/INS）
#    如你的 VCF 用 ALT=<DEL>/<INS> 但 INFO 里也有 SVTYPE，下面表达式即可工作。
#    若个别文件缺少 INFO/SVTYPE，可改成 -i 'ALT~"^<DEL>$"' 之类（见注释）。
FILTERED_VCF="${TMPDIR}/input.${SVT}.vcf.gz"
{ grep '^#' "$VCF_IN"; grep "${SVT}" "$VCF_IN"; } > "$FILTERED_VCF"

# 1) 从 TE 注释 TSV 生成 CHROM POS END TE_ID TE_CLASS
awk -F'\t' 'NR==1{print "CHROM\tPOS\tEND\tTE_ID\tTE_CLASS"; next} {split($1,a,"[:-]"); print a[2]"\t"a[3]"\t"a[3]"\t"$4"\t"$5}' \
  "$TE_TSV" > "$TMPDIR/te.table"

# 2) 去表头并用 faidx 顺序排序
tail -n +2 "$TMPDIR/te.table" > "$TMPDIR/te.table.nohdr"
bedtools sort -faidx "$FAI" -i "$TMPDIR/te.table.nohdr" > "$TMPDIR/te.table.sorted"

# 3) 生成用于 bcftools annotate 的 4 列（CHROM POS TE_ID TE_CLASS）
awk -v OFS='\t' '{print $1,$2,$4,$5}' "$TMPDIR/te.table.sorted" > "$TMPDIR/te.bed2"

# 4) 压缩索引
bgzip -c "$TMPDIR/te.bed2" > "$TMPDIR/te.bed2.gz"
tabix -s1 -b2 -e2 "$TMPDIR/te.bed2.gz"

# 5) 处理 header
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

# 6) 注释
OUT_ANN="${OUTPREFIX}.${SVT}.annotated.vcf.gz"
bcftools annotate \
  -a "$TMPDIR/te.bed2.gz" \
  -c CHROM,POS,INFO/TE_ID,INFO/TE_CLASS \
  -h "$HDR" \
  "$FILTERED_VCF" \
  -O z -o "$OUT_ANN"
tabix -p vcf "$OUT_ANN"

# 7) 仅保留带 TE_ID 的位点
OUT_TONLY="${OUTPREFIX}.${SVT}.TEonly.vcf.gz"
bcftools view -O z -o "$OUT_TONLY" -i 'INFO/TE_ID!=""' "$OUT_ANN"
tabix -p vcf "$OUT_TONLY"
rm "$OUT_ANN"
rm "$OUT_ANN.tbi"
echo "DONE."
echo "Annotated VCF : $OUT_ANN"
echo "TE-only VCF   : $OUT_TONLY"

# --- 备选：如果你的 VCF **没有** INFO/SVTYPE，而只用 ALT=<INS>/<DEL>:
# 将第 0 步替换为：
#   bcftools view -O z -o "$FILTERED_VCF" -i 'ALT~"^<'"$SVT"'>$"' "$VCF_IN"
#   tabix -p vcf "$FILTERED_VCF"


