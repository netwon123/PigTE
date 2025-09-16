awk -F'\t' '
  BEGIN{OFS="\t"}
  /^#/ {next}
  {
    ref=$4
    n=split($5, alts, ",")
    for(i=1;i<=n;i++){
      alt=alts[i]
      if (alt ~ /^</) continue
      if (length(ref) > length(alt)) {
        del=substr(ref, 2)                 
        gsub(/[^ACGTNacgtn]/, "N", del)   
        print ">SV-" $1 ":" $2
        print del
      }
    }
  }' panel_del.vcf > panel_del.fa

awk -F'\t' '
  BEGIN{OFS="\t"}
  /^#/ {next}
  {
    ref = $4
    n = split($5, alts, ",")
    for (i=1; i<=n; i++) {
      alt = alts[i]
      # 跳过符号性ALT（如 <INS>, <INS:ME> 等）
      if (alt ~ /^</) continue
      # 插入事件：ALT更长
      if (length(alt) > length(ref)) {
        ins = substr(alt, 2)                 # 去掉左锚定第1碱基
        gsub(/[^ACGTNacgtn]/, "N", ins)      # 可选：清理非常规字符
        print ">SV:" $1 "-" $2
        print ins
      }
    }
  }' panel_ins.vcf > panel_ins.fa

