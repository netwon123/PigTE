(
  head -n 1 susScr11.fa.out_filted08.bed
  awk  'NR>1{
    cov = $24         
    k   = $10             

    # 首次见到该key时记录顺序
    if (!(k in seen)) { seen[k] = ++n; keys[n] = k; bestcov[k] = -1 }

    # 选择规则：若已不是coverage=1，则遇到coverage=1就替换；否则选更大的coverage
    if (bestcov[k] != 1 && (cov == 1 || cov > bestcov[k])) { best[k] = $0; bestcov[k] = cov }
  }
  END{
    for (i=1; i<=n; i++) print best[keys[i]]
  }' OFS='\t' susScr11.fa.out_filted08.bed
) > TE_ref.dedup2.bed

awk '{print $5,$6,$7, $10"#"$11}' OFS='\t' TE_ref.dedup2.bed > TE_ref.dedup.bed
rm TE_ref.dedup2.bed

tail -n +2 TE_ref.dedup.bed > TE_ref.bed
rm TE_ref.dedup.bed
bedtools getfasta -fi ~/01_project_PGG/04_result/00_backone_genome/backone.fa -bed TE_ref.bed -name -fo TE_ref.tmp.fa
sed -E 's/^>([^:]+).*/>\1/' TE_ref.tmp.fa > TE_ref.fa
rm TE_ref.tmp.fa
