#!/usr/bin/bash
#SBATCH -p short -N 1 -n 2


# get the genome and transcriptome from FungiDB
mkdir -p genome
curl -o genome/JEL423_transcripts.fasta https://fungidb.org/common/downloads/release-48/BdendrobatidisJEL423/fasta/data/FungiDB-48_BdendrobatidisJEL423_AnnotatedTranscripts.fasta
curl -o genome/JEL423_genome.fasta https://fungidb.org/common/downloads/release-48/BdendrobatidisJEL423/fasta/data/FungiDB-48_BdendrobatidisJEL423_Genome.fasta
curl -o genome/JEL423_genes.gff https://fungidb.org/common/downloads/release-48/BdendrobatidisJEL423/gff/data/FungiDB-48_BdendrobatidisJEL423.gff

module load kallisto
kallisto index -i genome/JEL423_transcripts.idx genome/JEL423_transcripts.fasta -k 27

# if on the cluster - this is a copy from the cluster
mkdir -p input
rsync -a /bigdata/stajichlab/shared/projects/Chytrid/BdVirus/RNAseq_Bd_infected/snakemake_rnaseq_STAR/fastq/*.gz input/
