#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=2:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

### print start date and time
echo Job started on:
date -u

SCRIPTDIR=./Scripts/mQTL
InDir=./WGCNA/mQTL
OutDir=./WGCNA/mQTL/darkred.mQTL.results
covar_fact=Sex,TissueType,Plate,BraakStage
covar_num=Age,CellProportion
FilePrefix=darkred

mkdir -p $OutDir

module load R

for i in {1..22}
do
  echo "Running chr${i}..."
  plink --bfile ${InDir}/${FilePrefix} --recodeA --chr $i --out ${InDir}/${FilePrefix}_chr${i}
done

gcta64 --bfile ${InDir}/${FilePrefix} --make-grm-bin --out ${InDir}/${FilePrefix} --thread-num 16
gcta64 --grm ${InDir}/${FilePrefix} --pca --out ${InDir}/${FilePrefix}

Rscript ${SCRIPTDIR}/matrixQTL_eQTL_PrepareData.r $InDir $OutDir $covar_fact $covar_num $FilePrefix 


echo Job finished:
date -u
