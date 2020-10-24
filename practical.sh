#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -pe smp 8
#$ -N BIOS60132_Practical_Nine

#Load julia module
module load julia

#Set inputs
inPath="../Data/Practical_Nine"
gffFile="ref/GCF_000001405.39_GRCh38.p13_genomic.gff"
sampleList="SRR6408214 SRR6408215 SRR6408216"
outFile="geneCounts.csv"

#Run julia script to generate gene counts
julia ./countGenes.jl "$inPath" "$outFile" "$gffFile" "$sampleList"