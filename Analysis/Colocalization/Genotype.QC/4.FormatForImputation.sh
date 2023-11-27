## format files for use with Michegan Imputation Server

## EXECUTION
# sh SNPArray/preprocessing/4_formatForimputation.sh <population> <SNP ref file>
# where 
# <population > is 3 letter code for super population state ALL for no subsetting by population
# <SNP ref file> is an input file of 
# script needs to be executed from <git repo>/array/

## REQUIRES the following variables in config file
# PROCESSDIR, IMPUTEDIR, FILEPREFIX

## REQUIRES the following software
# plink, perl,

## INPUT
#  # binary plink files following prelim QC

## OUTPUT
# vcf files split by chr for upload to michegan imputation server

population=$1
refFile=$2

cd ${PROCESSDIR}

mkdir -p ImputationInput_${FILEPREFIX}

cd ImputationInput_${FILEPREFIX}

mkdir -p ${population}



## use tool to check data prior to upload https://www.well.ox.ac.uk/~wrayner/tools/
## for All use 1000G
cd ${population}

echo ${population}
## subset samples
#if [ $population != "ALL" ]
#then
    #plink --bfile ${PROCESSDIR}/QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd --keep ${PROCESSDIR}/merge1KG_${FILEPREFIX}/${population}Samples.txt --maf 0.05 --out ${FILEPREFIX}_QCd_${population} --make-bed
#else
	cp ${PROCESSDIR}/QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.bim ${FILEPREFIX}_QCd_${population}.bim
	cp ${PROCESSDIR}/QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.bed ${FILEPREFIX}_QCd_${population}.bed
	cp ${PROCESSDIR}/QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.fam ${FILEPREFIX}_QCd_${population}.fam
#fi

## liftover to hg19 for imputation
#plink --bfile ${FILEPREFIX}_ImputationInput/${FILEPREFIX}_QCd_${population} --update-map ${GSAREF}liftoverhg19.txt 3 --make-bed --out ${FILEPREFIX}_QCd_hg19 

plink --bfile ${FILEPREFIX}_QCd_${population} --freq --out ${FILEPREFIX}_QCd_${population}

#if [[ $(basename ${refFile}) == HRC* ]] ;
#then
#perl ${SCRIPTDIR}/HRC-1000G-check-bim.pl -b ${FILEPREFIX}_QCd_${population}.bim -f ${FILEPREFIX}_QCd_${population}.frq -r ${refFile} -g --hrc
#else 
perl ${SCRIPTDIR}/Pre-Imputation/HRC-1000G-check-bim.pl -b ${FILEPREFIX}_QCd_${population}.bim -f ${FILEPREFIX}_QCd_${population}.frq -r ${refFile} -g --1000g
#fi

#pwd
sed -i  '1 s+'"${FILEPREFIX}_${population}_QCd"'+'"${PROCESSDIR}/QCoutput_${FILEPREFIX}/${FILEPREFIX}_${population}_QCd"'+' Run-plink.sh
sed -i  '6,$ s+--make-bed+--recode vcf+' Run-plink.sh
sh Run-plink.sh
#for file in *.vcf; do awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${file} > with_chr_${file}; vcf-sort with_chr_${file} | bgzip -c > ${file}.gz;done
for file in *.vcf; do vcf-sort ${file} | bgzip -c > ${file}.gz;done
 
rm *.vcf
rm *.txt *.log
rm *.b* *.f*

