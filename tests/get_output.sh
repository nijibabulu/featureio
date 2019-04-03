#! /bin/bash

fasta_extract fasta/GCF_000744065.1_ASM74406v1_genomic.fna NZ_BBIY01000160.1 | fasta_substr /dev/stdin 20 10  | tail -n 1 > output/GCF_000744065.1_ASM74406v1_genomic.fna_NZ_BBIY01000160.1.20-30
fasta_stat fasta/GCF_000744065.1_ASM74406v1_genomic.fna  | head -n -2 > output/GCF_000744065.1_ASM74406v1_genomic.fna-lengths
fasta_stat fasta/random.fa | head -n -2 > output/random.fa-lengths
