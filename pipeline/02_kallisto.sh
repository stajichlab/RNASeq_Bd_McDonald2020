#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 8gb --out logs/kallisto.log


CPUS=$SLURM_CPUS_ON_NODE

if [ ! $CPUS ]; then
    CPUS=1
fi

module load kallisto
DB=JEL423_transcripts.idx
RESULTS=results
mkdir -p $RESULTS
tail -n +2 units.tsv | while read SAMPLE REP FQ1 FQ2
do
	kallisto quant -i $DB --single -l 300 --sd 40 --bias -t $CPU -o $RESULTS/$SAMPLE.$REP $FQ1
done
