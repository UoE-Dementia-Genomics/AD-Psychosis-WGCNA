#!/bin/sh
#!/bin/bash
#SBATCH -A Research_Project1 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --array=1-22 ##equal to chromosomes


### print start date and time
echo Job started on:
date -u 

###############################################################

## Set pvOutputThreshold > 0 and pvOutputThreshold.cis = 0 (or use Matrix_eQTL_engine) 
##   to perform eQTL analysis without using gene/SNP locations. Associations significant at the pvOutputThreshold level are be recorded in output_file_name and in the returned object.

## Set pvOutputThreshold = 0 and pvOutputThreshold.cis > 0 
##   to perform eQTL analysis for local gene-SNP pairs only. Local associations significant at pvOutputThreshold.cis level will be recorded in output_file_name.cis and in the returned object.

## Set pvOutputThreshold > 0 and pvOutputThreshold.cis > 0 
##   to perform eQTL analysis with separate p-value thresholds for local and distant eQTLs. Distant and local associations significant at corresponding thresholds are recorded in output_file_name and output_file_name.cis   respectively and in the returned object. 
##   In this case the false discovery rate is calculated separately for these two sets of eQTLs.
## for more information see ?Matrix_eQTL_main in R 

###############################################################

######

InDir=./darkgreen
SCRIPTDIR=./Analysis/mQTL
chr=$SLURM_ARRAY_TASK_ID
FilePrefix=darkgreen
Dist=1000000
cis_pval=1e-5            ## between 0 and 1
trans_pval=1e-5          ## between 0 and 1
trans_cross_chr=no    ## yes or no



Rscript ${SCRIPTDIR}/matrixQTL_eQTL_main.r $InDir $chr $FilePrefix $Dist $cis_pval $trans_pval "$trans_cross_chr"


echo Job finished:
date -u
