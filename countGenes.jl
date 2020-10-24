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
println("Importing packages...")
using BioAlignments
using GenomicFeatures
using DataStructures
using DataFrames 
using XAM
using GFF3
using CSV

#Retrieve inputs
outPath=ARGS[1]
gffPath=ARGS[2]
samplePath=ARGS[3]
gffFile=gffPath * "/" * ARGS[4]
bamList=Array{String}(undef, length(5:length(ARGS)))
baiList=Array{String}(undef, length(5:length(ARGS)))
for i in 5:length(ARGS)
    index=i-4
    bamList[index]=samplePath * "/" * ARGS[i] * "_sorted.bam"
    baiList[index]=samplePath * "/" * ARGS[i] * "_sorted.bam.bai"
end

#Count total mRNAs from gff
println("Counting mRNAs from gff...")
gffReader=open(GFF3.Reader, gffFile)
total=0
for record in gffReader
    if(GFF3.featuretype(record) == "mRNA")
        #Count total mRNAs
        global total=total+1
    end
end

#Prepare to count mRNAs from bam files
rows=total+1
cols=length(bamList)+1
counts=Array{String}(undef, rows, cols)
counts[1,1]="mRNA"
#Count features for each bam file
for sample in 1:length(bamList)
    #Add sample to matrix header
    colIndex=sample+1
    global counts[1,colIndex]=ARGS[colIndex-1]
    #Output sample to status message
    outMsg="Counting features for sample " * ARGS[sample] * "..."
    println(outMsg)
    #Retrieve input files
    bamReader=open(BAM.Reader, bamList[sample], index=baiList[sample])
    gffReader=open(GFF3.Reader, gffFile)
    #(Re)set row counter
    rowIndex=2
    #Loop over features
    for feature in gffReader
        if(GFF3.featuretype(feature) == "mRNA")
            #Count current feature
            i=0
            for record in eachoverlap(bamReader, feature)
                i=i+1
            end
            #Add features to matrix by first sample
            if(colIndex == 2)
                tag="mRNA" * rowIndex
                global counts[rowIndex,1]=tag
            end
            #Add counts to matrix
            global counts[rowIndex,colIndex]=string(i)
            #Update row index by feature
            rowIndex=rowIndex+1
        end
    end
end

#Write counts table to csv
println("Writing results to file...")
CSV.write(outPath, DataFrame(counts), header=false)