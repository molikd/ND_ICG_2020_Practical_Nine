#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -N BIOS60132_Practical_Nine

#Load julia module
module load julia

#Set inputs
outFile="geneCounts.csv"
gffPath="../Data/Practical_Nine/ref"
samplePath="../Data/Practical_Nine/results/bam"
gffFile="GCF_000001405.39_GRCh38.p13_genomic.gff"
sampleList="SRR6408214 SRR6408215 SRR6408216"

#Run julia script to generate gene counts
julia ./countGenes.jl "$outFile" "$gffPath" "$samplePath" "$gffFile" "$sampleList"