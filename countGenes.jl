#!usr/bin/env julia

#Script to calculate the sequence diversity of a fasta file
#Usage: julia ./countGenes.jl
#Usage Ex: julia ./countGenes.jl

#Import and add necessary packages
#import Pkg
#Pkg.add("BioAlignments")
#Pkg.add("GenomicFeatures")
#Pkg.add("DataStructures")
#Pkg.add("DataFrames")
#Pkg.add("XAM")
#Pkg.add("GFF3")
#Pkg.add("CSV")

#Import packages
using BioAlignments
using GenomicFeatures
using DataStructures
using DataFrames 
using XAM
using GFF3
using CSV

#Retrieve inputs
inPath=ARGS[1]
outPath=ARGS[2]
inputGFF=inPath * /ref/ * ARGS[3]
bamList=Array{String}(undef, length(4:length(ARGS)))
baiList=Array{String}(undef, length(4:length(ARGS)))
for i in 4:length(ARGS)
    index=i-3
    bamList[index]=inPath * "/results/bam/" * ARGS[i] * ".bam"
    baiList[index]=inPath * "/results/bam/" * ARGS[i] * ".bai"
end

#Retrieve bam and gff files
gffReader=open(GFF3.Reader, inputGFF)
for file in bamList
    bamReader=open(BAM.Reader, file)
end

#Count total mRNAs from gff
total=0
for record in gffReader
    if(GFF3.featuretype(record) == "mRNA")
        total=total+1
    end
end

#Prepare to count mRNAs from bam files
rows=length(total)+1
cols=length(bamList+1)
counts=zeros(rows,cols)
rowIndex=2
colIndex=2
counts[1,1]="mRNA"
#Count features for each bam file
for sample in 1:length(bamList)
    bamReader=open(BAM.Reader, bamList[sample], index=baiList[sample])
    for feature in gffReader
        if(GFF3.featuretype(feature) == "mRNA")
            i=0
            for record in eachoverlap(bamReader, feature)
                i=i+1
            end
            #Add feature and counts to matrix
            counts[rowIndex,1]=feature
            counts[rowIndex,colIndex]=i
            #Update row index by feature
            rowIndex=rowIndex+1
        end
    end
    #Add sample to matrix header
    counts[1,colIndex]=ARGS[colIndex+2]
    #Update column index by sample
    colIndex=colIndex+1
end

#Write counts table to csv
CSV.write(outPath, DataFrame(counts), writeheader=false)