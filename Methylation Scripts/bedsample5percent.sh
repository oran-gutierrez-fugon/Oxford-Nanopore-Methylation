#!/bin/bash

# Set the percentage of samples you want (0.5% in this case)
percentage_samples=0.5

# Input and output file paths
input_file="/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/NP4-3cpgUCSC_SORTED.bed"
output_file="/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/NP4-3cpgUCSC_SORTED_sampled.bed"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Input file not found."
  exit 1
fi

# Get the total number of lines in the input file
total_lines=$(wc -l < "$input_file")

# Calculate the number of samples based on the percentage
num_samples=$((total_lines * 5 / 1000)) # Since 0.5% is 5/1000

# Check if the number of lines is less than the number of samples requested
if [ "$total_lines" -lt "$num_samples" ]; then
  echo "Number of samples requested is greater than the total number of lines in the input file."
  exit 1
fi

# Use awk to sample the lines randomly
awk -v num_samples="$num_samples" -v total_lines="$total_lines" 'BEGIN { srand(); } { if (rand() < (num_samples/total_lines)) { print; num_samples--; } total_lines--; }' "$input_file" > "$output_file"

echo "Sampling completed. $num_samples random samples saved to $output_file"
