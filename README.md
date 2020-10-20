This is reanalysis of McDonald et 2020 - https://pubmed.ncbi.nlm.nih.gov/31892375/
using Kallisto and DEseq2

On the UCRHPCC you can run this
```
git clone https://github.com/stajichlab/RNASeq_Bd_McDonald2020
cd RNASeq_Bd_McDonald2020
mkdir fastq
ln -s /bigdata/stajichlab/shared/projects/Chytrid/BdVirus/RNAseq_Bd_infected/snakemake_rnaseq_STAR/fastq/*.gz fastq
sbatch pipeline/01_initialize_db.sh  
# wait for the above to finish running
# check with squeue -u $USER
sbatch pipeline/02_kallisto.sh
# wait for the above to finish running
# check with squeue -u $USER
Rscript Rscripts/kallisto_profile.R
```
