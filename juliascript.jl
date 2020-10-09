#!/usr/bin/env julia


#modules in use
println("Loading packages...")
using GenomicFeatures
using BioAlignments
using DataStructures
using XAM
using GFF3
using DataFrames
using Statistics
using CSV
using ArgParse


#Parse arguments
println("parsing arguments...")
s = ArgParseSettings()
@add_arg_table! s begin
    "--gff"
    "--sorted_bam_1"
    "--bam_index_1"
    "--sample_name_1"
    "--sorted_bam_2"
    "--bam_index_2"
    "--sample_name_2"
    "--sorted_bam_3"
    "--bam_index_3"
    "--sample_name_3"
    "--output"
end

parsed_args = parse_args(ARGS, s)


#set arguments as variables
gff_file = parsed_args["gff"]
sorted_bam_1 = parsed_args["sorted_bam_1"]
bam_index_1 = parsed_args["bam_index_1"]
sample_name_1 = parsed_args["sample_name_1"]
sorted_bam_2 = parsed_args["sorted_bam_2"]
bam_index_2 = parsed_args["bam_index_2"]
sample_name_2 = parsed_args["sample_name_2"]
sorted_bam_3 = parsed_args["sorted_bam_3"]
bam_index_3 = parsed_args["bam_index_3"]
sample_name_3 = parsed_args["sample_name_3"]
output_csv = parsed_args["output"]


#function to create a datafram with gene names and counts for one sample
println("creating function...")
function gene_count_df(gff_file::String, sorted_bam::String, bam_index::String, sample_name::String)

    #readers
    gff_reader = open(GFF3.Reader, gff_file)
    bam_reader = open(BAM.Reader, sorted_bam, index = bam_index)

    #=test to see if each overlap is working
    for feature in gff_reader
        for record in eachoverlap(bam_reader, feature)
        println(record)
        end
    end
    =#

    #=check to see if counting is working
    for feature in gff_reader
        if(GFF3.featuretype(feature) == "mRNA")
            println(feature)
            i=0 
            for record in eachoverlap(bam_reader, feature)
                i=i+1
            end
            println(i)
        end
    end
    =#

    #make DataFrame
    df1 = DataFrame(feature = GFF3.Record[], count = Int[])


    #fill dataframe with feature and count
    for feature in gff_reader
        if(GFF3.featuretype(feature) =="mRNA")
            i=0
            for record in eachoverlap(bam_reader, feature)
                i = i+1
            end
            push!(df1, (feature,i))
        end
    end

    #=create column of rnaIDs based on "ID"
    #create array of gene names
    rna_IDs = String[]
    for i in 1:nrow(df1)
        push!(rna_IDs, GFF3.attributes(df1[i,1], "ID")[1])
    end
    =#

    #create column of geneIDs based on "Parent"
    parent_genes = String[]
    for i in 1:nrow(df1)
        push!(parent_genes, GFF3.attributes(df1[i,1], "Parent")[1])
    end

    #add genes to dataframe as column
    df1.genes = parent_genes

    #remove column 1 that has all the extraneous information
    df1a = select(df1, Not(:feature))

    #groupby genes
    gdf1a = groupby(df1a, :genes)


    #average the counts for duplicate genes
    gene_count_df = combine(gdf1a, :count => mean)


    #name the columns
    column_names = ["Gene", sample_name]

    #convert column names to symbols
    column_names_symbols = Symbol[]
    for i in 1:length(column_names)
        symb = Symbol(column_names[i])
        push!(column_names_symbols, symb)
    end

    names!(gene_count_df, column_names_symbols)
    
    #return dataframe
    return gene_count_df
end

#create dataframe for each sample
println("evaluating first sample...")
Sample1_df = gene_count_df(gff_file, sorted_bam_1, bam_index_1, sample_name_1)

println("evaluating second sample...")
Sample2_df = gene_count_df(gff_file, sorted_bam_2, bam_index_2, sample_name_2)

println("evaluating third sample...")
Sample3_df = gene_count_df(gff_file, sorted_bam_3, bam_index_3, sample_name_3)


#join three dataframes on gene name
println("joining dataframes...")
Final_df = join(Sample1_df, Sample2_df, Sample3_df, on = :Gene, kind = :outer)

#print out dataframe (takes a long time)
#println(Final_df)

#save to csv file
CSV.write(output_csv, Final_df)
println("read count dataframe saved as csv")