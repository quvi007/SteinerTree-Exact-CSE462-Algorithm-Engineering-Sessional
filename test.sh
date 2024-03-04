#!/bin/bash

# Input and output file names
input_file="testcases/instance001.gr"
output_file="output.csv"

# Use awk to extract the required tokens and save them to the CSV file
awk 'NR==2 {printf "%s,", $2} NR==3 {printf "%s,", $2} /Terminals / {printf "%s,", $2}' "$input_file" > "$output_file"

# Add a newline character at the end of the output file
echo >> "$output_file"
