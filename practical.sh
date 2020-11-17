#!/bin/bash


module load julia 

julia GeneCount.jl --opt1 SRR12782308_sorted.bam --opt2 SRR12782308_sorted.bam.bai --opt3 SRR12782309_sorted.bam --opt4 SRR12782309_sorted.bam.bai --opt5 SRR12782310_sorted.bam --opt6 SRR12782310_sorted.bam.bai --opt7 Homo_sapiens.GRCh37.87.gff3
