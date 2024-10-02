#!/bin/bash

input_dir=$1
sample=$2
output_dir=$3

mkdir -p "$output_dir"

zcat "${input_dir}/${sample}_R1_bismark_bt2_pe.deduplicated.CX_report.txt.gz" | awk -v sample="$sample" -v outdir="$output_dir" '
{
    file = outdir "/" sample "_R1_bismark_bt2_pe.deduplicated." $6 "_report." $1 ".txt";
    print $0 >> file;
    files[file] = 1
}
END {
    for (file in files) {
        close(file);
        system("gzip " file);
    }
}'