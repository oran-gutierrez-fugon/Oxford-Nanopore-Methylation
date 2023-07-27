#!/bin/bash

reference_file="/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/reference.txt"
bed_file="NP4-3cpg.bed"
output_file="NP4-3cpg_CorrectChr.bed"

# Read the reference file and store the data in an associative array
declare -A reference_data
while read -r col1 col2; do
  reference_data["$col1"]="$col2"
done < "$reference_file"

# Process the bed file and create the output
while read -r col1 rest; do
  new_col1="${reference_data[$col1]}"
  if [[ -n $new_col1 ]]; then
    echo -e "$new_col1\t$rest"
  fi
done < "$bed_file" | cut -f 2- > "$output_file"

echo "Output file saved as: $output_file"
