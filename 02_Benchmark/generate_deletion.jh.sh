#!/usr/bin/env bash
# 生成每个 BAM 的独立 jh 脚本

# 配置
MELT_JAR=~/software/MELTv2.2.2/MELT.jar
BED=/public/home/baoqi/02_project_panTE_new/03_benchmark/01_library/del_TE.bed
REF=~/01_project_PGG/04_result/00_backone_genome/backone.fa
WORKDIR=$(pwd)

for bam in *.bam; do
    sample="${bam##*/}"
    sample="${sample%.bam}"
    jhfile="${sample}.jh"

    cat > "$jhfile" <<EOF
#!/bin/sh
#JSUB -J nucmer
#JSUB -n 8
#JSUB -q normal
#JSUB -o out.%J
#JSUB -e err.%J
#JSUB -R "span[hosts=1]"

source /public/jhinno/unischeduler/conf/unisched

# 创建样本工作目录
mkdir -p ${WORKDIR}/${sample}

# 1) MELT Deletion-Genotype
java -Xmx150G -jar ${MELT_JAR} Deletion-Genotype \\
    -w ${WORKDIR}/${sample} \\
    -bamfile ${WORKDIR}/${bam} \\
    -bed ${BED} \\
    -h ${REF}

# 2) 生成 mergelist
ls ${WORKDIR}/${sample}/*.tsv > ${WORKDIR}/${sample}/del_list.txt

# 3) MELT Deletion-Merge
java -Xmx150G -jar ${MELT_JAR} Deletion-Merge \\
    -bed ${BED} \\
    -mergelist ${WORKDIR}/${sample}/del_list.txt \\
    -h ${REF} \\
    -o ${WORKDIR}/${sample}

# 4) 重命名合并后的 VCF
mv ${WORKDIR}/${sample}/DEL.final_comp.vcf ${WORKDIR}/${sample}/${sample}.del.vcf

# 5) 清理中间 tsv
rm ${WORKDIR}/${sample}/*.tsv
EOF

    echo "Generated script: ${jhfile}"
done
