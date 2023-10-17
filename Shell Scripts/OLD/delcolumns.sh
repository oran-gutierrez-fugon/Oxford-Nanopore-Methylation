#!/bin/bash

# Input and output file paths
input_file="/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/testUCSC_SORTED_sampled.bed"
output_file="/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/testUCSC_SORTED_1st8columns.bed"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Input file not found."
  exit 1
fi

# Use awk to print only the first 8 columns and save to the output file
awk -v OFS="\t" '{NF=8} 1' "$input_file" > "$output_file"

echo "Columns after the 8th column removed. Result saved to $output_file"
