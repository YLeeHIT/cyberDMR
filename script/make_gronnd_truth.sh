#!/bin/bash

set -euo pipefail

if [ $# -lt 7  ]; then
    echo "Usage: $0 cyberDMR_result.bed metilene_result.bed BSmooth_result.bed HOME_result.bed MethyLasso_result.bed DiffMethylTools_result.bed OUTDIR [MIN_OVERLAP_BP]"
    echo "Example: $0 cyberDMR_result.bed metilene_result.bed BSmooth_result.bed HOME_result.bed MethyLasso_result.bed DiffMethylTools_result.bed result 1"
    exit 1
fi

A="$1"; B="$2"; C="$3"; D="$4"; E="$5"; F="$6"
OUTDIR="$7"
MINBP="${8:-1}"  # 最小重叠碱基数，默认 1bp

mkdir -p "$OUTDIR"/{sorted,pairwise,consensus,logs}

# 检查 bedtools
if ! command -v bedtools >/dev/null 2>&1; then
    echo "[ERROR] bedtools not found. Please install bedtools." >&2
    exit 2
fi

tool_name() {
    local f="$1"
    f="${f##*/}"              # 去掉路径
    f="${f%_result.bed}"      # 去掉后缀
    echo "$f"        
}

sort_bed() {
    local in="$1"
    local out="$2"
    awk 'BEGIN{OFS="\t"} NF>=3 && $2<$3 {print $1,$2,$3}' "$in" \
        | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "$out"
}

# 生成排序后的临时文件
S1="$OUTDIR/sorted/$(tool_name "$A").sorted.bed"
S2="$OUTDIR/sorted/$(tool_name "$B").sorted.bed"
S3="$OUTDIR/sorted/$(tool_name "$C").sorted.bed"
S4="$OUTDIR/sorted/$(tool_name "$D").sorted.bed"
S5="$OUTDIR/sorted/$(tool_name "$E").sorted.bed"
S6="$OUTDIR/sorted/$(tool_name "$F").sorted.bed"

echo "[INFO] Sorting inputs ..."
sort_bed "$A" "$S1"
sort_bed "$B" "$S2"
sort_bed "$C" "$S3"
sort_bed "$D" "$S4"
sort_bed "$E" "$S5"
sort_bed "$F" "$S6"

# 文件和工具名数组
FILES=("$S1" "$S2" "$S3" "$S4" "$S5" "$S6")
NAMES=($(tool_name "$A") $(tool_name "$B") $(tool_name "$C") $(tool_name "$D") $(tool_name "$E") $(tool_name "$F"))

# 两两求交集
echo "[INFO] Generating pairwise overlaps ..."
for ((i=0;i<6;i++)); do
    for ((j=i+1;j<6;j++)); do
        n1="${NAMES[$i]}"
        n2="${NAMES[$j]}"
        out="$OUTDIR/pairwise/${n1}_vs_${n2}.overlap.bed"
        bedtools intersect -a "${FILES[$i]}" -b "${FILES[$j]}" -wao \
            | awk -v OFS="\t" -v MINBP="$MINBP" '
                $NF >= MINBP {
                    chrA=$1; sA=$2; eA=$3; chrB=$4; sB=$5; eB=$6;
                    if (chrA==chrB) {
                        s = (sA > sB ? sA : sB);
                        e = (eA < eB ? eA : eB);
                        if (e > s) print chrA, s, e;                                                                 
                    }
                                                
                }' \
            | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
            | bedtools merge -i - > "$out"
        echo "  -> $out"
    done
done

echo "[INFO] Building consensus (>=2 tools) ..."
bedtools multiinter -i "$S1" "$S2" "$S3" "$S4" "$S5" "$S6" > "$OUTDIR/consensus/multiinter.raw.bed"

awk 'BEGIN{OFS="\t"} NR>1 && $4>=2 {print $1,$2,$3}' "$OUTDIR/consensus/multiinter.raw.bed" \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
    | bedtools merge -i - > "$OUTDIR/consensus/consensus_ge2.bed"

echo "  -> $OUTDIR/consensus/consensus_ge2.bed"

# 汇总统计
echo "[INFO] Summarizing ..."
{
    echo -e "Pair\tRegions\tTotalBp"
    for f in "$OUTDIR"/pairwise/*.overlap.bed; do
        regions=$(wc -l < "$f")
        bp=$(awk '{sum+=($3-$2)} END{print sum+0}' "$f")
        echo -e "$(basename "$f")\t${regions}\t${bp}"
    done
} > "$OUTDIR/logs/pairwise_summary.tsv"
