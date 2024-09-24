#!/bin/bash

prefix=$1
output_dir=$2

mkdir -p "$output_dir"

zcat ${prefix}_R1_bismark_bt2_pe.deduplicated.CX_report.txt.gz | awk -v pfx="$prefix" -v outdir="$output_dir" '
{
    file = outdir "/" pfx "_R1_bismark_bt2_pe.deduplicated." $6 "_report." $1 ".txt";
    print $0 >> file;
    files[file] = 1
}
END {
    for (file in files) {
        close(file);
        system("gzip " file);
    }
}'