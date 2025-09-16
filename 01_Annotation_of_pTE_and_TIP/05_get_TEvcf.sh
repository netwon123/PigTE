awk -F'\t' 'NR==1{print "CHROM\tPOS\tEND\tTE_ID\tTE_CLASS"; next}  \
{                       \
  split($1,a,"[:-]"); \
  print a[2]"\t"a[3]"\t"a[3]"\t"$4"\t"$5 \
}' Best_TE_annotation_per_SV_ref-TE.tsv  > Best_TE_annotation_per_SV_ref-TE.tsv.table 


tail -n +2 Best_TE_annotation_per_SV_ref-TE.tsv.table > Best_TE_annotation_per_SV_ref-TE.tsv.table.bed

 bedtools sort -faidx ~/01_project_PGG/04_result/00_backone_genome/backone.fa.fai \
 -i Best_TE_annotation_per_SV_ref-TE.tsv.table.bed \
 > Best_TE_annotation_per_SV_ref-TE.tsv.table.sort.bed
 
 
 
 awk '{print $1,$2,$4,$5}' OFS='\t' Best_TE_annotation_per_SV_ref-TE.tsv.table.sort.bed > Best_TE_annotation_per_SV_ref-TE.tsv.table.bed2

 bgzip Best_TE_annotation_per_SV_ref-TE.tsv.table.bed2
tabix -s1 -b2 -e2 Best_TE_annotation_per_SV_ref-TE.tsv.table.bed2.gz

bcftools  annotate -a Best_TE_annotation_per_SV_ref-TE.tsv.table.bed2.gz  -c CHROM,POS,INFO/TE_ID,INFO/TE_CLASS -h ../header.txt ../panel_del.vcf -o panel_del.annotated.vcf -O v

awk '/^#/ {print; next} $8 ~ /TE_ID=/' panel_del.annotated.vcf > panel_del.TEonly.vcf