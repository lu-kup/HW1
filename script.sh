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