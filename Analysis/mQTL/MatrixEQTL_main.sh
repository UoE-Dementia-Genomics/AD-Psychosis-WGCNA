#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --array=1-22


### print start date and time
echo Job started on:
date -u

InDir=./WGCNA/mQTL/darkred.mQTL.results
chr=$SLURM_ARRAY_TASK_ID
FilePrefix=darkred
Dist=1e6
cis_pval=0.05
trans_pval=0.05
SCRIPTDIR=./Scripts/mQTL

module load R

Rscript ${SCRIPTDIR}/matrixQTL_eQTL_main.r $InDir $chr $FilePrefix $Dist $cis_pval $trans_pval 


echo Job finished:
date -u
