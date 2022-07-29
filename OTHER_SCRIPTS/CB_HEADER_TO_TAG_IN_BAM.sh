#!/usr/bin/env bash

# Load modules
module load sambamba

# Extract header
sambamba view -H Sample.bam > tmpHeader.txt

# Add CB to tag
sambamba view -t 23 Sample.bam | awk '/^S/ {N=split($1,n,"_")-1;print $0 "\tCB:Z:" n[N]}' > tmp.sam

# Create sam file with header
cat tmpHeader.txt tmp.sam > tmp_with_header.sam
rm tmpHeader.sam tmp.sam

# Create bam file with CB tag
samtools view -@ 23 -bS tmp_with_header.sam > Sample_With_CB.bam
samtools index Sample_With_CB.bam




