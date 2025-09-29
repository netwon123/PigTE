#!/usr/bin/env bash
# generate_melt_jh.sh
# 遍历当前目录下所有 BAM 文件，生成独立 MELT jh 作业脚本

MELT_JAR=~/software/MELTv2.2.2/MELT.jar
REF=~/01_project_PGG/04_result/00_backone_genome/backone.fa
ZIP_DIR=/public/home/baoqi/02_project_panTE_new/03_benchmark/01_library/melt_zip/zips
GENE_BED=/public/home/baoqi/02_project_panTE_new/03_benchmark/01_library/genes_melt_format.bed
WORKDIR=$(pwd)

for bam in *.bam; do
    sample="${bam##*/}"
    sample="${sample%.bam}"
    jhfile="${sample}.jh"
    outdir="${WORKDIR}/${sample}_melt"

    cat > "$jhfile" <<EOF
#!/bin/sh
#JSUB -J ${sample}_melt
#JSUB -n 8
#JSUB -q normal
#JSUB -o out.%J
#JSUB -e err.%J
#JSUB -R "span[hosts=1]"

source /public/jhinno/unischeduler/conf/unisched

mkdir -p ${outdir}

# 1) Preprocess
java -jar ${MELT_JAR} Preprocess -bamfile ${WORKDIR}/${bam} -h ${REF}

# 2) IndivAnalysis
for zip in ${ZIP_DIR}/*zip; do
    java -jar ${MELT_JAR} IndivAnalysis \
        -bamfile ${WORKDIR}/${bam} \
        -w ${outdir} \
        -h ${REF} \
        -t \$zip \
        -c 50 \
        -r 150
done

# 3) GroupAnalysis
for zip in ${ZIP_DIR}/*zip; do
    java -jar ${MELT_JAR} GroupAnalysis \
        -discoverydir ${outdir} \
        -h ${REF} \
        -t \$zip \
        -r 150 \
        -n ${GENE_BED} \
        -w ${outdir}
done

# 4) Genotype
for zip in ${ZIP_DIR}/*zip; do
    java -jar ${MELT_JAR} Genotype \
        -bamfile ${WORKDIR}/${bam} \
        -w ${outdir} \
        -p ${outdir} \
        -h ${REF} \
        -t \$zip
done

# 5) MakeVCF
for zip in ${ZIP_DIR}/*zip; do
    java -jar ${MELT_JAR} MakeVCF \
        -genotypingdir ${outdir} \
        -p ${outdir} \
        -o ${outdir} \
        -w ${outdir} \
        -t \$zip \
        -h ${REF}
done

# 6) 压缩和索引 VCF
cd ${outdir}
for vcf in *.vcf; do
    bgzip -@ 8 \$vcf
done

for gz in *.gz; do
    tabix -f \$gz
done

ls *.vcf.gz > merge.list
bcftools merge -l merge.list -Oz -o ${sample}.ins.vcf.gz --threads 8
mv ${sample}.ins.vcf.gz ../
cd ../
rm -r ${outdir}
EOF

    echo "Generated MELT job script: ${jhfile}"
done

