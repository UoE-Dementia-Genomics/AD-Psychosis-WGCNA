
population=$1
refFile=$2

cd ${PROCESSDIR}

mkdir -p ImputationInput_${FILEPREFIX}

cd ImputationInput_${FILEPREFIX}

mkdir -p ${population}

## for All use 1000G
cd ${population}

echo ${population}

cp ${PROCESSDIR}/QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.bim ${FILEPREFIX}_QCd_${population}.bim
cp ${PROCESSDIR}/QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.bed ${FILEPREFIX}_QCd_${population}.bed
cp ${PROCESSDIR}/QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.fam ${FILEPREFIX}_QCd_${population}.fam

## liftover to hg19 for imputation
#plink --bfile ${FILEPREFIX}_ImputationInput/${FILEPREFIX}_QCd_${population} --update-map ${GSAREF}liftoverhg19.txt 3 --make-bed --out ${FILEPREFIX}_QCd_hg19 

plink --bfile ${FILEPREFIX}_QCd_${population} --freq --out ${FILEPREFIX}_QCd_${population}

perl ${SCRIPTDIR}/HRC-1000G-check-bim.pl -b ${FILEPREFIX}_QCd_${population}.bim -f ${FILEPREFIX}_QCd_${population}.frq -r ${refFile} -g --1000g

sed -i  '1 s+'"${FILEPREFIX}_${population}_QCd"'+'"${PROCESSDIR}/QCoutput_${FILEPREFIX}/${FILEPREFIX}_${population}_QCd"'+' Run-plink.sh
sed -i  '6,$ s+--make-bed+--recode vcf+' Run-plink.sh
sh Run-plink.sh

for file in *.vcf; do vcf-sort ${file} | bgzip -c > ${file}.gz;done
 
rm *.vcf
rm *.txt *.log
rm *.b* *.f*

