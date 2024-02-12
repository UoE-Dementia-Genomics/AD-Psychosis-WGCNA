#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=1:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err


ScriptDir=./EWCE/Scripts

dir_10x=./ROSMAP/SC
pheno_file=./ROSMAP/SC/filtered_column_metadata.csv
TargetGenes_file=./darkred.methyl-expr.nominally.sig.genes.1kb.csv
rep_=10000
n_cores=16
generate_ctd=T
ctd_file=./ctd_scROS.rda
out_prefix=./EWCE/darkred

Rscript ${ScriptDir}/EWCE.CelltypeEnrichment.R $dir_10x $pheno_file $TargetGenes_file $rep_ $n_cores $generate_ctd $ctd_file $out_prefix

