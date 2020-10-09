#!/bin/bash
#$ -M sweaver4@nd.edu
#$ -m abe
#$ -N sweaver4_Practical_Nine

module load julia

julia juliascript.jl --gff "genome.gff" \
	--sorted_bam_1 "sample_1/SRR6408410.sorted.bam" \
	--bam_index_1 "sample_1/SRR6408410.sorted.bam.bai" \
    --sample_name_1 "SRR6408410" \
    --sorted_bam_2 "sample_2/SRR6408239.sorted.bam" \
	--bam_index_2 "sample_2/SRR6408239.sorted.bam.bai" \
    --sample_name_2 "SRR6408239" \
    --sorted_bam_3 "sample_3/SRR6408314.sorted.bam" \
	--bam_index_3 "sample_3/SRR6408314.sorted.bam.bai" \
    --sample_name_3 "SRR6408314" \
    --output "Gene_counts_3_samples.csv"