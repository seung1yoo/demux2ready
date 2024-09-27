#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: $0 <Control> <Case> <Output_Dir> <Context>"
    exit 1
fi

control=$1
case=$2
output_dir=$3
context=$4

mkdir -p "${output_dir}/${context}"

pushd "${output_dir}/${context}" || exit

metilene_input.pl --in1 ../../${control}_R1_bismark_bt2_pe.deduplicated.bedGraph.${context}.gz \
                  --in2 ../../${case}_R1_bismark_bt2_pe.deduplicated.bedGraph.${context}.gz \
                  --h1 "$control" --h2 "$case"

metilene -t 4 -a "$control" -b "$case" metilene_${control}_${case}.input | \
    sort -V -k1,1 -k2,2n > metilene_${control}_${case}.output

metilene_output.pl -q metilene_${control}_${case}.output \
                   -o metilene_${control}_${case}.output.filter \
                   -a "$control" -b "$case"

bedtools intersect -a metilene_${control}_${case}.output.filter_qval.0.05.out \
                   -b ../../anno/MANE.GRCh38.v1.3.summary.chr.sorted.bed \
                   -wa -loj > metilene_${control}_${case}.output.filter_qval.0.05.out.anno

sed -i "1ichr\tstart\tstop\tq-value\tmean_methylation_delta\tcpg_count\tmean_${control}\tmean_${case}\tanno_chr\tanno_start\tanno_stop\tanno_gene\tanno_score\tanno_strand" \
    metilene_${control}_${case}.output.filter_qval.0.05.out.anno

# 8. bedtools intersect로 전체 주석 추가
bedtools intersect -a metilene_${control}_${case}.output \
                   -b ../../anno/MANE.GRCh38.v1.3.summary.chr.sorted.bed \
                   -wa -loj > metilene_${control}_${case}.output.anno

sed -i "1ichr\tstart\tstop\tq-value\tmean_methylation_delta\tcpg_count\tp-MWU\tp-2D_KS\tmean_${control}\tmean_${case}\tanno_chr\tanno_start\tanno_stop\tanno_gene\tanno_score\tanno_strand" \
    metilene_${control}_${case}.output.anno

popd || exit

echo "Processing completed successfully for $control and $case in ${output_dir}/${context}"
