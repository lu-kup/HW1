
hisat2 -q -p 6 -x references/indexed -1 reads_trimmed/SRR11647705_1_val_1.fq.gz -2 reads_trimmed/SRR11647705_2_val_2.fq.gz -S SRR11647705.sam --verbose
hisat2 -q -p 6 -x references/indexed -1 reads_trimmed/SRR11647706_1_val_1.fq.gz -2 reads_trimmed/SRR11647706_2_val_2.fq.gz -S SRR11647706.sam --verbose
hisat2 -q -p 6 -x references/indexed -1 reads_trimmed/SRR11647707_1_val_1.fq.gz -2 reads_trimmed/SRR11647707_2_val_2.fq.gz -S SRR11647707.sam --verbose
hisat2 -q -p 6 -x references/indexed -1 reads_trimmed/SRR11647708_1_val_1.fq.gz -2 reads_trimmed/SRR11647708_2_val_2.fq.gz -S SRR11647708.sam --verbose
hisat2 -q -p 6 -x references/indexed -1 reads_trimmed/SRR11647709_1_val_1.fq.gz -2 reads_trimmed/SRR11647709_2_val_2.fq.gz -S SRR11647709.sam --verbose
hisat2 -q -p 6 -x references/indexed -1 reads_trimmed/SRR11647710_1_val_1.fq.gz -2 reads_trimmed/SRR11647710_2_val_2.fq.gz -S SRR11647710.sam --verbose
mkdir repeated_outputs
mv *.sam repeated_outputs