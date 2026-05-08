#!/bin/bash
prefetch SRR22316123
prefetch SRR22316122
prefetch SRR22316118
prefetch SRR22316117

fastq-dump --split-3 SRR22316123.sra
fastq-dump --split-3 SRR22316122.sra
fastq-dump --split-3 SRR22316118.sra
fastq-dump --split-3 SRR22316117.sra

ls ref_data/GSE218147/*_1.fastq >1
ls ref_data/GSE218147/*_2.fastq >2
paste 1 2 >data_tcf.config
rm 1 2


#### clean data  ####
# trim------
mkdir ./trim
mkdir ./trim_qc
cat data_tcf.config | while read id 
do 
arr=($id)
fq1=${arr[0]}
fq2=${arr[1]}
trim_galore --phred33 --gzip --no_report_file -q 20 --trim-n -o ref_data/GSE218147/trim --paired $fq1 $fq2
done

#### mapping_mm10####
mkdir ./mapping_mm10
mkdir ./mapping_mm10/bowtie2_summary

cat data_tcf.config | while read id
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
-1 ref_data/GSE218147/trim/${name1}'_val_1.fq.gz' \
-2 ref_data/GSE218147/trim/${name2}'_val_2.fq.gz' \
-S ref_data/GSE218147/mapping_mm10/${outname}_bowtie2.sam \
&> ref_data/GSE218147/mapping_mm10/bowtie2_summary/${outname}_bowtie2.txt

done



#### filter and transformat ####
mkdir ./mapping_mm10/bam

cat data_tcf.config | while read id 
do 
arr=($id)
fq1=${arr[0]}
fq2=${arr[1]}
##echo $fq1,$fq2

basename=$(basename "$fq1")
name=${basename%%_*}
echo $name

## Filter and keep the mapped read pairs
samtools view -bS -F 0x04 ref_data/GSE218147/mapping_mm10/${name}_bowtie2.sam -o ref_data/GSE218147/mapping_mm10/bam/${name}_bowtie2.mapped.bam

done


# mkdir ./peakCalling_mm10
macs3 callpeak -t ref_data/GSE218147/mapping_mm10/bam/SRR22316129_bowtie2.mapped.bam ref_data/GSE218147/mapping_mm10/bam/SRR22316130_bowtie2.mapped.bam \
-c ref_data/GSE218147/mapping_mm10/bam/SRR22316124_bowtie2.mapped.bam ref_data/GSE218147/mapping_mm10/bam/SRR22316125_bowtie2.mapped.bam  ref_data/GSE218147/mapping_mm10/bam/SRR22316126_bowtie2.mapped.bam \
-g mm --bdg -f BAMPE -n Tcf7_macs3 --outdir ref_data/GSE218147/peakCalling_mm10 -q 0.05 2>ref_data/GSE218147/peakCalling_mm10/Tcf7_Peak_summary.txt

# bedGraph to bigWig
bedGraphToBigWig  peakCalling_mm10/Tcf7_macs3_treat_pileup.bdg  mm10.fai  Tcf7.bw


