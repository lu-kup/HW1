for f in mapping_outputs/*.sorted.bam; do
    clipping_profile.py -s "PE" -i "$f" -o outputs/${f%.sorted.bam}/clipping_profile
done

for f in mapping_outputs/*.sorted.bam; do
    junction_annotation.py -r raw_data/hg38_GENCODE_V47.bed -i "$f" -o outputs/${f%.sorted.bam}/junction_annotation
done

for f in mapping_outputs/*.sorted.bam; do
    infer_experiment.py -r raw_data/hg38_GENCODE_V47.bed -i "$f" > outputs/${f%.sorted.bam}/infer_experiment.out
done

for f in mapping_outputs/*.sorted.bam; do
    inner_distance.py -r raw_data/hg38_GENCODE_V47.bed -i "$f" -o outputs/${f%.sorted.bam}/inner_distance
done

for f in mapping_outputs/*.sorted.bam; do
    geneBody_coverage.py -r raw_data/hg38_GENCODE_V47.bed -i "$f" -o outputs/${f%.sorted.bam}/gene_body_coverage
done