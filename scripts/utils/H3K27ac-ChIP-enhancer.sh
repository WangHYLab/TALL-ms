#!/bin/bash

## analysis for H3K27ac, then compare with Bcl11a ChIP-seq data to find Bcl11a-bound enhancers

# phred check----
ls -l SRR31588264.fastq| head -1 |awk '{print $9}' >tmp.txt
cat tmp.txt | while read id
do
echo $id
# zcat $id | head -100 >tmp2.fastq
cat $id | head -100 >tmp2.fastq
sh rna-analysis/fq_qual_type.sh tmp2.fastq
done
rm ./tmp.txt
rm ./tmp2.fastq

#### clean data  ####
## trim------
mkdir ./trim
trim_galore --phred33 --gzip --no_report_file -q 20 --trim-n -o ref_data/GSE283392/trim SRR31588264.fastq


#### mapping mm10####
mkdir ./mapping
mkdir ./mapping/bowtie2_summary

bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 50 \
	-x /home/public_data/bowtie2_ref/mm10/mm10 \
	-U ref_data/GSE283392/trim/SRR31588264_trimmed.fq.gz \
	-S ref_data/GSE283392/mapping/SRR31588264_bowtie2.sam \
	&> ref_data/GSE283392/mapping/bowtie2_summary/SRR31588264_bowtie2.txt



#### filter and transformat ####
mkdir ./mapping/bam

## Filter and keep the mapped read pairs
samtools view -bS -F 0x04 ref_data/GSE283392/mapping/SRR31588264_bowtie2.sam \
	-o ref_data/GSE283392/mapping/bam/SRR31588264_bowtie2.mapped.bam

## Convert into bed file format
bedtools bamtobed -i ref_data/GSE283392/mapping/bam/SRR31588264_bowtie2.mapped.bam -bedpe > ref_data/GSE283392/mapping/bed/SRR31588264_bowtie2.bed


mkdir ./peakCalling/
macs3 predictd -i ref_data/GSE283392/mapping/bam/SRR31588264_bowtie2.mapped.bam \
	-f BAM --outdir ref_data/GSE283392/peakCalling -g mm
# out is 209

macs3 callpeak -t ref_data/GSE283392/mapping/bam/SRR31588264_bowtie2.mapped.bam \
	-g mm --bdg -f BAM -n DN2_H3K27ac \
	 --nomodel --broad --broad-cutoff 0.1 --extsize 209  --shift -104 \
	--outdir ref_data/GSE283392/peakCalling -q 0.05 \
	2>ref_data/GSE283392/peakCalling/DN2_H3K27ac_Peak_summary.txt

## enhancer
cd ref_data/GSE283392/peakCalling/
# 5.1: extract TSS points for all transcripts
GTF="/home/public_data/mm10/genes.gtf"
awk 'BEGIN {OFS="\t"; FS="\t"} 
  $3 == "transcript" {
    # extract basic info
    chrom = $1;
    start = $4;
    end = $5;
    strand = $7;
    
    # recognize gene_id, gene_name, transcript_id
    split($9, attr, ";");
    gene_id = "";
    gene_name = "";
    transcript_id = "";
    
    for (i in attr) {
      gsub(/^[ \t]+|[ \t]+$/, "", attr[i]);  
      if (attr[i] ~ /^gene_id/) {
        split(attr[i], tmp, "\"");
        gene_id = tmp[2];
      }
      if (attr[i] ~ /^gene_name/) {
        split(attr[i], tmp, "\"");
        gene_name = tmp[2];
      }
      if (attr[i] ~ /^transcript_id/) {
        split(attr[i], tmp, "\"");
        transcript_id = tmp[2];
      }
    }
    
    # make TSS point
    tss = (strand == "+") ? start : end;
    
    # generate unique ID for each TSS point
    id = (gene_name != "") ? gene_name "_" transcript_id : gene_id "_" transcript_id;
    
    # output TSS point
    print chrom, tss, tss + 1, id, ".", strand;
}' $GTF | sort -k1,1 -k2,2n > all_tss_points.bed


# 5.2: create TSS regions (±2.5kb)
# fasta
samtools faidx /home/public_data/bowtie2_ref/mm10/genome.fa
cut -f1,2 /home/public_data/bowtie2_ref/mm10/genome.fa.fai > chrom.sizes

GENOME_SIZE="mm"
CHROM_SIZES='chrom.sizes'
bedtools slop -i all_tss_points.bed -g chrom.sizes -b 2500 > tss_regions_2.5kb.bed

# 5.3: filter H3K27ac wide peak (FDR<0.1)
awk '$9 <= 0.1' DN2_H3K27ac_peaks.broadPeak > H3K27ac_fdr0.1.bed


# 5.4: recognize candidate enhancers (exclude TSS regions)
bedtools subtract -a DN2_H3K27ac_peaks.broadPeak -b tss_regions_2.5kb.bed -A > H3K27ac_enhancers_candidate.bed

# 5.5: filter small enhancers (optional)
awk '{if ($3-$2 >= 200) print}' H3K27ac_enhancers_candidate.bed > H3K27ac_enhancers_filtered.bed

# 5.6: recognize Bcl11a-bound enhancers
bedtools intersect -a cuttag-analysis/peakCalling_mm10/B_merge_macs3_peaks.narrowPeak \
  -b H3K27ac_enhancers_filtered.bed -wa -u > Bcl11a_bound_enhancers.bed

# 5.7: add annotation
awk 'BEGIN{OFS="\t"; count=1} 
  {print $1,$2,$3,"Enhancer_"count++,$5,$6,$7,$8,$9}' Bcl11a_bound_enhancers.bed > Bcl11a_bound_enhancers_annotated.bed