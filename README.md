# ND_ICG_2020_Practical_Nine
Practical Nine for Introduction to Computational Genomics

**note that I can't put the bam and gff files in the repo because they are too big. But I downloaded fasta files using the SRA toolkit, aligned them using bwa, converted to bam, sorted, etc, and then used with the practical.sh and juliascript.jl files that are here**

This practical is intended to quiz you

Download this repository, and create a Julia script that takes three samples of this [RNAseq data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108383) and the [human genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39) and creates and prints a three column array of the counts of genes per sample. You are allowed to use the aligner of your choice to align trimmed fatq files to the genome file, but you must use Julia to make the table of samples x genes. Create a practical.sh script which calls your Julia script with the included fasta file.
