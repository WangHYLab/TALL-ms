#!/bin/bash

# ENSMUST00000109514.7
EXON_LENGTH=5284
# EXON_REGION="chr11:24163145-24168429" #real
EXON_REGION="11:24162845-24168729" #flank

OUTPUT_FILE="WES-analysis/WES.bcl11a_exon4_analysis.csv"
echo "Sample_ID,Reads_Count,Total_Reads,TPM,RPKM" > $OUTPUT_FILE

cd WES-analysis/out/bam

while IFS= read sample
do
	arr=($sample)
	fq1=${arr[0]}
	fq2=${arr[1]}
	# echo $fq1
	name=${fq1#*WES/}
	sample=${name%%_*}
	echo $sample

    # 
    BAM_FILE="${sample}_bqsr.bam"

    # 
    READS_COUNT=$(samtools view -c "$BAM_FILE" "$EXON_REGION")

    # 
    TOTAL_READS=$(samtools view -c "$BAM_FILE")

    # TPM
    if [ "$TOTAL_READS" -ne 0 ]; then
        TPM=$(echo "scale=4; ($READS_COUNT / ($EXON_LENGTH / 1000)) / ($TOTAL_READS / 1000000)" | bc)
    else
        TPM=0
    fi

    # RPKM
    if [ "$TOTAL_READS" -ne 0 ]; then
        RPKM=$(echo "scale=4; ($READS_COUNT / ($EXON_LENGTH / 1000)) / ($TOTAL_READS / 1000000)" | bc)
    else
        RPKM=0
    fi

    echo "$sample,$READS_COUNT,$TOTAL_READS,$TPM,$RPKM" >> $OUTPUT_FILE

done < WES-analysis/WES_data.config

echo "Analysis complete. Results saved to $OUTPUT_FILE."