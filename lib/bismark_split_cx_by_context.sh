#!/bin/bash
# best usage : bismark_split_cx_by_context.sh \
#                  {input_dir -> demux2ready/bismark} \
#                  {sample -> A7_1} \
#                  {output_dir -> demux2ready/bismark/A7_1_split_cx_by_context}
#                  {is_dedup -> deduplicated.}

input_dir=$1
sample=$2
output_dir=$3
is_dedup=$4

mkdir -p "$output_dir"

zcat "${input_dir}/${sample}_R1_bismark_bt2_pe.${is_dedup}CX_report.txt.gz" | awk -v sample="$sample" -v outdir="$output_dir" -v is_dedup="$is_dedup" '
{
    file = outdir "/" sample "_R1_bismark_bt2_pe." is_dedup $6 "_report." $1 ".txt";
    print $0 >> file;
    files[file] = 1
}
END {
    for (file in files) {
        close(file);
        system("gzip " file);
    }
}'
