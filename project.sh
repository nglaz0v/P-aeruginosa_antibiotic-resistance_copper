#!/bin/bash -x

# [SRX9370194: GSM4862160: PAO1_MIC_Exp5_t=72 (DNA-seq); Pseudomonas aeruginosa; OTHER](https://www.ncbi.nlm.nih.gov/sra/SRR12905311)
# Study: Do Copper enhance antibiotic resistance (AR) in resistant bacteria and/or induce AR in sensitive bacteria? [DNA-seq]
# Organism: [Pseudomonas aeruginosa (Синегнойная палочка)](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=287)
# Reference genome: [ASM676v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006765.1)
# GEO Accession: [GSM4862160](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4862160)
# SRA: [SRP288731](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP288731)
# BioProject: [PRJNA672287](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA672287)
# GEO: [GSE160188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160188)

ORGANISM=PAeruginosa
GENOME=GCF_000006765.1
RUN=SRR12905311

ABRICATE=~/Apps/abricate/bin/abricate
QUALIMAP=~/Apps/qualimap_v2.3/qualimap


RESULT_DIR=/mnt/c/Users/Asus/BioGen
[ -d ${RESULT_DIR} ] && rm -rf ${RESULT_DIR} && mkdir ${RESULT_DIR}

#mkdir -p project && cd project

mkdir -p reads ref

# Скачивание ридов в формате FASTQ
[ -f reads/${RUN}_1.fastq.gz ] || wget -P reads ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/011/SRR12905311/SRR12905311_1.fastq.gz
[ -f reads/${RUN}_2.fastq.gz ] || wget -P reads ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR129/011/SRR12905311/SRR12905311_2.fastq.gz

# Скачивание референсного генома для Pseudomonas aeruginosa в формате FASTA и GFF3
[ -f ${GENOME}.zip ] || curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000006765.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000006765.1.zip" -H "Accept: application/zip"
unzip ${GENOME}.zip
mv ncbi_dataset/data/${GENOME}/${GENOME}_ASM*_genomic.fna ref/${ORGANISM}.fasta
mv ncbi_dataset/data/${GENOME}/genomic.gff ref/${ORGANISM}.gff
rm -rf ncbi_dataset README.md

rm -rf quality_check && mkdir -p quality_check
# Тримминг
mkdir -p trimmed
fastp --overrepresentation_analysis --trim_poly_g \
    -f 12 -t 2 \
    -h quality_check/${RUN}.fastp.html -j quality_check/${RUN}.fastp.json --thread `nproc` \
    -i reads/${RUN}_1.fastq.gz -I reads/${RUN}_2.fastq.gz -o trimmed/${RUN}_1.fastq.gz -O trimmed/${RUN}_2.fastq.gz
firefox quality_check/${RUN}.fastp.html

FASTQ_DIR=trimmed  # reads || trimmed
# Оценка качества ридов
fastqc -o quality_check $FASTQ_DIR/*.fastq.gz ../${RUN}.fastq.gz
firefox quality_check/${RUN}_*_fastqc.html
rm -rf multiqc_data multiqc_report.html
multiqc quality_check $FASTQ_DIR -o quality_check
firefox quality_check/multiqc_report.html || mv quality_check ${RESULT_DIR}

#mkdir -p assembly
# Creating a genome assembly
#spades.py -o assembly/spades-original --careful -1 $FASTQ_DIR/${RUN}_1.fastq.gz -2 $FASTQ_DIR/${RUN}_2.fastq.gz --threads `nproc`
# Assembly quality assessment
#quast.py -o assembly/quast assembly/spades-original/scaffolds.fasta --threads `nproc`
#firefox assembly/quast/report.pdf || mv assembly/quast/report.pdf ${RESULT_DIR}/assembly_report.pdf
#${ABRICATE} assembly/spades-original/scaffolds.fasta --threads `nproc` > abricate_report.tsv

# Выравнивание ридов на геном
mkdir -p index mappings
# generate index for mapping
bwa index ref/${ORGANISM}.fasta
#bowtie2-build ref/${ORGANISM}.fasta index/${ORGANISM}
# mapping reads
bwa mem ref/${ORGANISM}.fasta $FASTQ_DIR/${RUN}_1.fastq.gz $FASTQ_DIR/${RUN}_2.fastq.gz > mappings/${ORGANISM}.sam -t `nproc`
#bowtie2 -x index/${ORGANISM} -1 $FASTQ_DIR/${RUN}_1.fastq.gz -2 $FASTQ_DIR/${RUN}_2.fastq.gz -S mappings/${ORGANISM}.sam --threads `nproc`

# converting output from SAM to BAM
samtools view -b -o mappings/${ORGANISM}.bam mappings/${ORGANISM}.sam
rm mappings/${ORGANISM}.sam
# show statistics
samtools flagstat mappings/${ORGANISM}.bam
# sorting by BAM coordinate (good practice)
samtools sort -@ 8 -o mappings/${ORGANISM}.sorted.bam mappings/${ORGANISM}.bam

# building index for BAM and FASTA (to access the files faster)
samtools index mappings/${ORGANISM}.sorted.bam
samtools faidx ref/${ORGANISM}.fasta

# collect BAM statistics and generate plots
rm -rf bam_stats
samtools stats mappings/${ORGANISM}.sorted.bam > mappings/${ORGANISM}.sorted.bam.bc
plot-bamstats -p bam_stats/ mappings/${ORGANISM}.sorted.bam.bc
firefox bam_stats/index.html
${QUALIMAP} bamqc -bam mappings/${ORGANISM}.sorted.bam
firefox mappings/${ORGANISM}.sorted_stats/qualimapReport.html

# SNP calling
mkdir -p variants
#bcftools mpileup -f ref/${ORGANISM}.fasta mappings/${ORGANISM}.sorted.bam | bcftools call -mv -o calls.vcf && cat calls.vcf | grep -v ^## | head -5 && bcftools view -i '%QUAL>=20' calls.vcf
freebayes -p 1 -f ref/${ORGANISM}.fasta mappings/${ORGANISM}.sorted.bam > variants/${ORGANISM}.freebayes.vcf  # '-p 1' - specifies the ploidy level: E.Coli are haploid

# collect VCF statistics and generate plots
rm -rf vcf_stats
bgzip -f variants/${ORGANISM}.freebayes.vcf
tabix -p vcf variants/${ORGANISM}.freebayes.vcf.gz
bcftools stats -F ref/${ORGANISM}.fasta -s - variants/${ORGANISM}.freebayes.vcf.gz > variants/${ORGANISM}.freebayes.vcf.gz.stats
plot-vcfstats -p vcf_stats/ variants/${ORGANISM}.freebayes.vcf.gz.stats
firefox vcf_stats/summary.pdf
#bgzip -d variants/${ORGANISM}.freebayes.vcf.gz

# Отображение результатов
igv -g ref/${ORGANISM}.fasta mappings/${ORGANISM}.sorted.bam variants/${ORGANISM}.freebayes.vcf.gz ref/${ORGANISM}.gff

# Удаление ненужных данных
mv mappings/${ORGANISM}.sorted_stats ./qualimap_report
[ -d ${RESULT_DIR} ] && mv bam_stats vcf_stats qualimap_report ${RESULT_DIR}
mkdir -p data && mv ref/${ORGANISM}.fasta mappings/${ORGANISM}.sorted.bam* variants/${ORGANISM}.freebayes.vcf.gz GSE160188_core.vcf.gz ref/${ORGANISM}.gff data
rm -rf variants mappings index assembly trimmed ref
