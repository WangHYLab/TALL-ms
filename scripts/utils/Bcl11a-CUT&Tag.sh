#!/bin/bash

ls rawdata/cuttag/*_1*fq.gz >1
ls rawdata/cuttag/*_2*fq.gz >2
paste 1 2 >data.config
rm 1 2

#### qc ####
mkdir cuttag-analysis/qc/
fastqc rawdata/cuttag/*.gz -o cuttag-analysis/qc/
multiqc cuttag-analysis/qc/. -o cuttag-analysis/qc/ -n multiqc_cuttag

#### phred check ####
#ls -l rawdata/cuttag/B1_L1_1.fq.gz| head -1 |awk '{print $9}' >tmp.txt
#cat tmp.txt | while read id
#do
#echo $id
#zcat $id | head -100 >tmp2.fastq
##cat $id | head -100 >tmp2.fastq
#sh /home/zhengjie/scripts/fq_qual_type.sh tmp2.fastq
#done
#rm ./tmp.txt
#rm ./tmp2.fastq

#### clean data  ####
## trim
mkdir cuttag-analysis/trim
mkdir cuttag-analysis/trim_qc
cat data.config | while read id 
do 
arr=($id)
fq1=${arr[0]}
fq2=${arr[1]}
trim_galore --phred33 --gzip --no_report_file -q 20 --trim-n -o cuttag-analysis/trim --paired $fq1 $fq2
done
fastqc cuttag-analysis/trim/*.gz -o cuttag-analysis/trim_qc/
multiqc cuttag-analysis/trim_qc/. -o cuttag-analysis/trim_qc/ -n multiqc_trim_cuttag

#### mapping ####
mkdir cuttag-analysis/mapping_mm10
mkdir cuttag-analysis/mapping_mm10/bowtie2_summary

cat data.config | while read id
do
arr=($id)
fq1=${arr[0]}
fq2=${arr[1]}
#echo $fq1,$fq2

basename=$(basename "$fq1")
outname=${basename%%_*}

name1=$(basename "$fq1")
name1=${name1%%.*}
name2=$(basename "$fq2")
name2=${name2%%.*}

## mouse
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 50 \
-x /home/public_data/bowtie2_ref/mm10/mm10 \
-1 cuttag-analysis/trim/${name1}'_val_1.fq.gz' \
-2 cuttag-analysis/trim/${name2}'_val_2.fq.gz' \
-S cuttag-analysis/mapping_mm10/${outname}_bowtie2.sam \
&> cuttag-analysis/mapping_mm10/bowtie2_summary/${outname}_bowtie2.txt

## E.coli
# bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --no-overlap --no-dovetail --phred33 -I 10 -X 700 -p 50 \
# -x /home/public_data/bowtie2_ref/E.coli/E.coli \
# -1 cuttag-analysis/trim/${name1}'_val_1.fq.gz' \
# -2 cuttag-analysis/trim/${name2}'_val_2.fq.gz' \
# -S cuttag-analysis/mapping_mm10/${outname}_bowtie2_E.coli.sam \
# &> cuttag-analysis/mapping_mm10/bowtie2_summary/${outname}_bowtie2_E.coli.txt

done




#### filter and transformat ####
mkdir cuttag-analysis/mapping_mm10/fragmentLen
mkdir cuttag-analysis/mapping_mm10/bed
mkdir cuttag-analysis/mapping_mm10/bam

cat data.config | while read id 
do 
arr=($id)
fq1=${arr[0]}
fq2=${arr[1]}
##echo $fq1,$fq2

basename=$(basename "$fq1")
name=${basename%%_*}
echo $name

## Filter and keep the mapped read pairs
samtools view -bS -F 0x04 cuttag-analysis/mapping_mm10/${name}_bowtie2.sam -o cuttag-analysis/mapping_mm10/bam/${name}_bowtie2.mapped.bam

## Convert into bed file format
bedtools bamtobed -i cuttag-analysis/mapping_mm10/bam/${name}_bowtie2.mapped.bam -bedpe > cuttag-analysis/mapping_mm10/bed/${name}_bowtie2.bed

## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
awk '$1==$4 && $6-$2 < 1000 {print $0}' cuttag-analysis/mapping_mm10/bed/${name}_bowtie2.bed > cuttag-analysis/mapping_mm10/bed/${name}_bowtie2.clean.bed

## Only extract the fragment related columns
cut -f 1,2,6 cuttag-analysis/mapping_mm10/bed/${name}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  > cuttag-analysis/mapping_mm10/bed/${name}_bowtie2.fragments.bed

## fragment size information ##
## Extract the 9th column from the alignment sam file which is the fragment length
samtools view -F 0x04 cuttag-analysis/mapping_mm10/${name}_bowtie2.sam |\
awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' |\
sort |\
uniq -c |\
awk -v OFS="\t" '{print $2, $1/2}' > cuttag-analysis/mapping_mm10/fragmentLen/${name}_fragmentLen.txt

## Sample Repeatability ####
binLen=500
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' cuttag-analysis/mapping_mm10/bed/${name}_bowtie2.fragments.bed |\
	sort -k1,1V -k2,2n |\
	uniq -c |\
	awk -v OFS="\t" '{print $2, $3, $1}' |\
	sort -k1,1V -k2,2n  > cuttag-analysis/mapping_mm10/bed/${name}_bowtie2.fragmentsCount.bin$binLen.bed


done


#### peak calling ####
mkdir ./peakCalling_mm10

## peak calling for each replicate
macs3 callpeak -t cuttag-analysis/mapping_mm10/bam/B1_bowtie2.mapped.bam \
-c cuttag-analysis/mapping_mm10/bam/IgG_bowtie2.mapped.bam \
-g mm --bdg -f BAMPE -n B1_macs3 --outdir cuttag-analysis/peakCalling_mm10 -q 0.05 2>cuttag-analysis/peakCalling_mm10/B1_Peak_summary.txt

macs3 callpeak -t cuttag-analysis/mapping_mm10/bam/B2_bowtie2.mapped.bam \
-c cuttag-analysis/mapping_mm10/bam/IgG_bowtie2.mapped.bam \
-g mm --bdg -f BAMPE -n B2_macs3 --outdir cuttag-analysis/peakCalling_mm10 -q 0.05 2>cuttag-analysis/peakCalling_mm10/B2_Peak_summary.txt

macs3 callpeak -t cuttag-analysis/mapping_mm10/bam/B3_bowtie2.mapped.bam \
-c cuttag-analysis/mapping_mm10/bam/IgG_bowtie2.mapped.bam \
-g mm --bdg -f BAMPE -n B3_macs3 --outdir cuttag-analysis/peakCalling_mm10 -q 0.05 2>cuttag-analysis/peakCalling_mm10/B3_Peak_summary.txt

## merge peaks from three replicates
macs3 callpeak -t cuttag-analysis/mapping_mm10/bam/B1_bowtie2.mapped.bam cuttag-analysis/mapping_mm10/bam/B2_bowtie2.mapped.bam cuttag-analysis/mapping_mm10/bam/B3_bowtie2.mapped.bam \
-c cuttag-analysis/mapping_mm10/bam/IgG_bowtie2.mapped.bam \
-g mm --bdg -f BAMPE -n B_merge_macs3 --outdir cuttag-analysis/peakCalling_mm10 -q 0.05 2>cuttag-analysis/peakCalling_mm10/Merge_Peak_summary.txt
