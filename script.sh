fastqc raw_data/*.fastq.gz -t 8 -o qc
trim_galore -o reads_trimmed -j 8 --quality 20 --length 20 --paired raw_data/SRR11647705_1.fastq.gz raw_data/SRR11647705_2.fastq.gz
trim_galore -o reads_trimmed -j 8 --quality 20 --length 20 --paired raw_data/SRR11647706_1.fastq.gz raw_data/SRR11647706_2.fastq.gz
trim_galore -o reads_trimmed -j 8 --quality 20 --length 20 --paired raw_data/SRR11647707_1.fastq.gz raw_data/SRR11647707_2.fastq.gz
trim_galore -o reads_trimmed -j 8 --quality 20 --length 20 --paired raw_data/SRR11647708_1.fastq.gz raw_data/SRR11647708_2.fastq.gz
trim_galore -o reads_trimmed -j 8 --quality 20 --length 20 --paired raw_data/SRR11647709_1.fastq.gz raw_data/SRR11647709_2.fastq.gz
trim_galore -o reads_trimmed -j 8 --quality 20 --length 20 --paired raw_data/SRR11647710_1.fastq.gz raw_data/SRR11647710_2.fastq.gz
mkdir qc_trimmed
fastqc reads_trimmed/*.fq.gz -t 8 -o qc_trimmed
gunzip references/GRCh38.primary_assembly.genome.fa.gz
hisat2-build references/GRCh38.primary_assembly.genome.fa references/indexed

hisat2 -q -p 8 -x references/indexed -1 reads_trimmed/SRR11647705_1_val_1.fq.gz -2 reads_trimmed/SRR11647705_2_val_2.fq.gz -S SRR11647705.sam --verbose
hisat2 -q -p 8 -x references/indexed -1 reads_trimmed/SRR11647706_1_val_1.fq.gz -2 reads_trimmed/SRR11647706_2_val_2.fq.gz -S SRR11647706.sam --verbose
hisat2 -q -p 8 -x references/indexed -1 reads_trimmed/SRR11647707_1_val_1.fq.gz -2 reads_trimmed/SRR11647707_2_val_2.fq.gz -S SRR11647707.sam --verbose
hisat2 -q -p 8 -x references/indexed -1 reads_trimmed/SRR11647708_1_val_1.fq.gz -2 reads_trimmed/SRR11647708_2_val_2.fq.gz -S SRR11647708.sam --verbose
hisat2 -q -p 8 -x references/indexed -1 reads_trimmed/SRR11647709_1_val_1.fq.gz -2 reads_trimmed/SRR11647709_2_val_2.fq.gz -S SRR11647709.sam --verbose
hisat2 -q -p 8 -x references/indexed -1 reads_trimmed/SRR11647710_1_val_1.fq.gz -2 reads_trimmed/SRR11647710_2_val_2.fq.gz -S SRR11647710.sam --verbose
mkdir mapping_outputs
mv *.sam mapping_outputs
samtools view -bS mapping_outputs/*.sam > mapping_outputs/*.bam

for f in mapping_outputs/*.sam; do
  samtools view -bS "$f" > "${f%.sam}.bam"
done

for f in mapping_outputs/*.bam; do
  samtools sort "$f" -o "${f%.bam}.sorted.bam"
done

samtools index -@ 8 -M mapping_outputs/*.sorted.bam

mkdir mapped_quality_assesment

for f in mapping_outputs/*.sorted.bam; do
  samtools coverage "$f" > "${f%.bam}_coverage.out"
done

mkdir outputs
mv mapping_outputs/*.out outputs

wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE_V47.bed.gz/download
cp download raw_data/hg38_GENCODE_V47.bed.gz
rm download

gunzip raw_data/hg38_GENCODE_V47.bed.gz

geneBody_coverage.py -r raw_data/hg38_GENCODE_V47.bed -i mapping_outputs/*.sorted.bam -o outputs/gene_body_coverage
inner_distance.py -r raw_data/hg38_GENCODE_V47.bed -i mapping_outputs/*.sorted.bam -o outputs/inner_distance
clipping_profile.py -s "PE" -i mapping_outputs/*.sorted.bam -o outputs/clipping_profile
junction_annotation.py -r raw_data/hg38_GENCODE_V47.bed -i mapping_outputs/*.sorted.bam -o outputs/junction_annotation
infer_experiment.py -r raw_data/hg38_GENCODE_V47.bed -i mapping_outputs/*.sorted.bam > outputs/infer_experiment.out

multiBamSummary bins --bamfiles mapping_outputs/*.sorted.bam -o bins_results.npz
plotCorrelation -in bins_results.npz -c spearman -p heatmap -o outputs/correlation_heatmap.png
plotPCA -in bins_results.npz -o outputs/pca.png

multiqc outputs

featureCounts -T 8 -a references/gencode.v49.annotation.gtf.gz -o outputs/feature_counts.out \
-s 0 -p -t exon -g gene_id mapping_outputs/*.sorted.bam

zcat references/gencode.v49.annotation.gtf.gz \
| awk '$3=="gene" { 
    match($0, /gene_name "([^"]+)"/, a); 
    match($0, /gene_type "([^"]+)"/, b); 
    match($0, /gene_id "([^"]+)"/, c); 
    print a[1] "," b[1] "," c[1]
}' > outputs/gene_table.out