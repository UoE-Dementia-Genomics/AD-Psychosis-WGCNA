#!/bin/bash


cd ${PROCESSDIR}

mkdir -p QCoutput_${FILEPREFIX}

##convert ped/map to bed
if [ ! -f ${RAWDATADIR}/${FILEPREFIX}.bim ] || [ ! -f ${RAWDATADIR}/${FILEPREFIX}.bed ] || [ ! -f ${RAWDATADIR}/${FILEPREFIX}.fam ]; then
	plink --file ${RAWDATADIR}/${FILEPREFIX} --maf 0.05 --make-bed --out  ${RAWDATADIR}/${FILEPREFIX}
fi

## some samples were run twice (two different brain regions) identify these:
king -b ${RAWDATADIR}/${FILEPREFIX}.bed --duplicate --prefix QCoutput_${FILEPREFIX}/${FILEPREFIX}_king

if [ -f QCoutput_${FILEPREFIX}/${FILEPREFIX}_king.con ]
then
	cut -f 1,2 QCoutput_${FILEPREFIX}/${FILEPREFIX}_king.con > QCoutput_${FILEPREFIX}/${FILEPREFIX}_duplicateSamples.tmp
	cut -f 3,4 QCoutput_${FILEPREFIX}/${FILEPREFIX}_king.con >> QCoutput_${FILEPREFIX}/${FILEPREFIX}_duplicateSamples.tmp
	sort QCoutput_${FILEPREFIX}/${FILEPREFIX}_duplicateSamples.tmp | uniq > QCoutput_${FILEPREFIX}/${FILEPREFIX}_duplicateSamples.txt
	rm QCoutput_${FILEPREFIX}/${FILEPREFIX}_duplicateSamples.tmp
fi

##if any duplicate samples calc missingness rate to identify those to exclude
if [ -s QCoutput_${FILEPREFIX}/${FILEPREFIX}_duplicateSamples.txt ]
then 
	plink --bfile ${RAWDATADIR}/${FILEPREFIX} --missing --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_duplicateSamples

	## use python script to identify duplicated with greatest missingness
	python ExcludeDuplicates.py QCoutput_${FILEPREFIX}/${FILEPREFIX}_king.con QCoutput_${FILEPREFIX}/${FILEPREFIX}_duplicateSamples.imiss QCoutput_${FILEPREFIX}/${FILEPREFIX}_dupsToExclude.txt

	## remove duplicates
	plink --bfile ${RAWDATADIR}/${FILEPREFIX} --remove QCoutput_${FILEPREFIX}/${FILEPREFIX}_dupsToExclude.txt --make-bed --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_1
else
	cp ${RAWDATADIR}/${FILEPREFIX}.bed QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_1.bed
	cp ${RAWDATADIR}/${FILEPREFIX}.bim QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_1.bim
	cp ${RAWDATADIR}/${FILEPREFIX}.fam QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_1.fam
fi	


## remove variants at the same position (i.e. triallelic)
awk '{if ($1 != 0) print $1":"$4}' QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_1.bim > pos.tmp
sort pos.tmp | uniq -d > dupLocs.txt
awk -F ":" '{print $1,$2-1,$2,"set1", "set2"}' dupLocs.txt > QCoutput_${FILEPREFIX}/${FILEPREFIX}_positionsExclude.txt

plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_1 --exclude range QCoutput_${FILEPREFIX}/${FILEPREFIX}_positionsExclude.txt --make-bed --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_2

rm pos.tmp
rm dupLocs.txt

## perform sex check on samples with enough data
plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_2 --mind 0.02 --check-sex --out QCoutput_${FILEPREFIX}/${FILEPREFIX}

## we keep samples which do not have sex predicted but check their predicted sex with other softwares such as GenomeStudio
## exclude mismatched samples
## If a samples doesn't have real sex, we can use the predicted sex for it
awk '{if ( $4 != $3 && $3 != 0 && $4 != 0 ) print $1,$2,$3,$4,$5,$6}'  QCoutput_${FILEPREFIX}/${FILEPREFIX}.sexcheck >> QCoutput_${FILEPREFIX}/${FILEPREFIX}_sexErrors.txt
plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_2 --remove QCoutput_${FILEPREFIX}/${FILEPREFIX}_sexErrors.txt --make-bed --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_3

## check for runs of homozygosity
awk '{if ($1 >= 1 && $1 <= 22) print $2}' QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_3.bim > QCoutput_${FILEPREFIX}/autosomalVariants.txt
plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_3 --extract QCoutput_${FILEPREFIX}/autosomalVariants.txt --maf 0.01 --hwe 0.00001 --mind 0.02 --geno 0.05 --indep-pairwise 5000 1000 0.2 --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_ld.auto

plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_3 --extract QCoutput_${FILEPREFIX}/${FILEPREFIX}_ld.auto.prune.in --het --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_roh
## exclude anyone with |Fhet| > 0.2
awk '{if ($6 > 0.2 || $6 < -0.2) print $1,$2}' QCoutput_${FILEPREFIX}/${FILEPREFIX}_roh.het > QCoutput_${FILEPREFIX}/${FILEPREFIX}_excessHet.txt
plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_3 --remove QCoutput_${FILEPREFIX}/${FILEPREFIX}_excessHet.txt --make-bed --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_4
#rm autosomalVariants.txt

##remove variants with 3+ alleles
plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_4 --biallelic-only strict  --make-bed --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_5


## filter sample and variant missingness, HWE, rare variants and exclude variants with no position
awk '{if ($1 == 0) print $2}' QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_5.bim > QCoutput_${FILEPREFIX}/${FILEPREFIX}_noLocPos.tmp
plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_5 --exclude QCoutput_${FILEPREFIX}/${FILEPREFIX}_noLocPos.tmp --maf 0.05 --hwe 0.00001 --mind 0.05 --geno 0.01 --make-bed --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd

### 9816 variants removed due to missing genotype data (--geno).
### --hwe: 10 variants removed due to Hardy-Weinberg exact test.
### 464 variants removed due to minor allele threshold(s)
### 284719 variants and 222 people pass filters and QC.

## write list of samples that passed QC for CNV calling
#cut -f 1,2 --delimiter=" " QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.fam > QCoutput_${FILEPREFIX}/${FILEPREFIX}_ID_Map.txt
cut -f 2 --delimiter=" " QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.fam > QCoutput_${FILEPREFIX}/${FILEPREFIX}_PastQCSamples.txt

## clean up intermediate files but keep log files
rm QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_*.b*
rm QCoutput_${FILEPREFIX}/${FILEPREFIX}_update_*.fam
 

## calc PCS within sample only
# LD prune
plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd --indep 50 5 1.5 --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.ld
plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd --extract QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.ld.prune.in --make-bed --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.ld.prune

mkdir -p GCTA_${FILEPREFIX}

gcta64 --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.ld.prune --make-grm-bin --autosome --out GCTA_${FILEPREFIX}/${FILEPREFIX}_QCd_GCTA
gcta64 --grm GCTA_${FILEPREFIX}/${FILEPREFIX}_QCd_GCTA --pca --out GCTA_${FILEPREFIX}/${FILEPREFIX}_QCd.pca

rm QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.ld.prune*

## extract SNP probes for comparison with DNAm data
#plink --bfile QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd --extract ${RefDir}/RSprobes.txt --recodeA --out QCoutput_${FILEPREFIX}/${FILEPREFIX}_DNAmSNPs


echo "General QC was done"
