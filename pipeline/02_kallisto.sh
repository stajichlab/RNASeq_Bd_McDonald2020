#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 8gb --out logs/kallisto.log

CPU=$SLURM_CPUS_ON_NODE

if [ ! $CPUS ]; then
    CPU=1
fi

module load kallisto
DB=genome/JEL423_transcripts.idx
RESULTS=results
mkdir -p $RESULTS
tail -n +2 units.tsv | while read SAMPLE REP FQ1 FQ2
do
# this could be sped up by doing this with parallel instead of serially but you get the results either way
	kallisto quant -i $DB --single -l 300 --sd 30 --bias -t $CPU -o $RESULTS/$SAMPLE.$REP $FQ1
done
