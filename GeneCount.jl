#!/usr/bin/env julia
import Pkg
Pkg.add("XAM")
using XAM, BioAlignments
using GFF3, GenomicFeatures
Pkg.add("ArgParse")
using DataFrames
using CSV
using ArgParse
using GenomeGraphs
function parse_commandline()


s = ArgParseSettings()
@add_arg_table s begin
        "--opt1"
                help = "BAM1 Reads"
        "--opt2"

                help = "BAM1 indicies"

        "--opt3"
                help = "BAM2 Reads"
	"--opt4"
                help = "BAM2 indicies"
        "--opt5"

                help = "BAM3 Reads"

        "--opt6"
                help = "BAM3 indicies"
	"--opt7"
		help = "GFF3 file"
end

return parse_args(s)
end

Reads = parse_commandline()

reader1 = open(BAM.Reader, READS["--opt1"], index = READS["--opt2"])
reader2 = open(BAM.Reader, READS["--opt3"], index = READS["--opt4"])
reader3 = open(BAM.Reader, READS["--opt6"], index = READS["--opt6"])


features = open(collect, GFF3.Reader, READS["--opt7"])

#keep mRNA features

filter!(x -> GFF3.featuretype(x) == "mRNA", features)

println("read gff3 file")


i = 0
for feature in features
        global i = i+1
end

MATRIX = zeros[3,i]
i = 0
for feature in features
	global i = i+1
	for record in eachoverlap(reader1,feature)
		global MATRIX(1,i) = MATRIX(1,i)+1
	end 
end 

i = 0                                                                            for feature in features
        global i = i+1
	for record in eachoverlap(reader2,feature)
		global MATRIX(2,i) = MATRIX(2,i)+1
        end
end		


i = 0                                                                            for feature in features
        global i = i+1
	for record in eachoverlap(reader3,feature)
		global MATRIX(3,i) = MATRIX(3,i)+1
	end
end

genecount = DataFrame(MATRIX)

CSV.write("genecount.csv",genecount,header = false)









	

