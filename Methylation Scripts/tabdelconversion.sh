#!/bin/bash

# Replace spaces with tabs in the input file and save the output to a new file
input_file="/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/NP4-3cpgUCSC_SORTED_sampled.bed"
output_file="/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/test2/NP4-3cpgUCSC_SORTED_tabdel.bed"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Input file not found."
  exit 1
fi

# Convert spaces to tabs using awk
awk 'BEGIN { OFS="\t"; } {gsub(/ +/, "\t"); print}' "$input_file" > "$output_file"

echo "Conversion completed. Output saved to $output_file"
