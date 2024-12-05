#!/bin/bash
# best usage : bismark_cx_to_bedCov.sh \
#                  {input_dir -> demux2ready/bismark/A7_1_split_cx_by_context} \
#                  {sample -> A7_1} \
#                  {context -> CG} \
#                  {min_depth -> 5}
#                  {is_dedup -> "deduplicated." or ""}

input_dir=$1
sample=$2
context=$3
min_depth=$4
is_dedup=$5

output_file="${input_dir}/${sample}_R1_bismark_bt2_pe.${is_dedup}bedCov.${context}.minDepth${min_depth}"

if [ -f "$output_file" ]; then
    rm "$output_file"
fi

if [ -f "${output_file}.gz" ]; then
    rm "${output_file}.gz"
fi

if [ -f "${output_file}.temp" ]; then
    rm "${output_file}.temp"
fi

find "$input_dir" -name "${sample}_*.${context}_report.*.txt.gz" | sort | while read -r file; do
    echo "Processing file: $file"

    zcat $file | awk '
    {
        depth = $4 + $5
        if (depth < '$min_depth') {
            next
        }
        
        new_col1 = $1
        new_col2 = $2
        new_col3 = $2
        if (depth > 0) {
            ratio = ($4 / depth) * 100
        } else {
            ratio = 0
        }
        meth_depth = $4 
        unmeth_depth = $5

        printf "%s\t%d\t%d\t%.2f\t%d\t%d\n", new_col1, new_col2, new_col3, ratio, meth_depth, unmeth_depth
    }' >> "${output_file}.temp"
done


echo "Sorting : ${output_file}.temp"
bedtools sort -i "${output_file}.temp" > "${output_file}"
rm "${output_file}.temp"

echo "Compressing : ${output_file}"
gzip "${output_file}"
echo "BedCov files combined and compressed into: ${output_file}.gz"
