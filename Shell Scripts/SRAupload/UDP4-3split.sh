#!/bin/bash

# Name of the original gzipped file
GZ_FILE="UDP4-3cat.fast5.gz"
# Prefix for the split files
SPLIT_PREFIX="UDP4-3cat_part_"
# Size of each split file (90 GB here to ensure it's less than 100 GB)
SPLIT_SIZE="90G"

# Step 1: Decompress the gzipped file
echo "Decompressing $GZ_FILE..."
gunzip $GZ_FILE

# Extract the original file name (without .gz)
ORIGINAL_FILE="${GZ_FILE%.gz}"

# Step 2: Split the file into smaller chunks
echo "Splitting $ORIGINAL_FILE into smaller chunks..."
split --bytes=$SPLIT_SIZE $ORIGINAL_FILE $SPLIT_PREFIX

# Step 3: Compress the split files
echo "Compressing the split files..."
for file in ${SPLIT_PREFIX}*; do
    echo "Compressing $file..."
    gzip $file
done

echo "Splitting and compression completed."
