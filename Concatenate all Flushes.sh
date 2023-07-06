find $BAM_DIR -name '*.bam' | {
    read firstbam
    samtools view -h "$firstbam"
    while read bam; do
        samtools view "$bam"
    done
} > UDP4-3concat.sam

mv "UDP4-3concat.sam" "UDP4-3concat.bam"

# Print the message after the code finishes running
echo "It is Done"