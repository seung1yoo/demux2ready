#bedtools map -a target_bait.bed -b ../bismark/Donor_split_cx_by_context/Donor_R1_bismark_bt2_pe.deduplicated.bedCov.CG.gz -c 4,5,6 -o count,sum,sum | awk 'BEGIN {FS="\t"; OFS="\t"} {depth = $5 + $6; if (depth > 0) {ratio = $5 / depth * 100} else {ratio = "NA"}; print $1, $2, $3, $4, $5, $6, ratio}' | sed '1ichr\tstart\tend\tcytosine_count\tmethyl_sum\tunmethyl_sum\tmethyl_level' | head  > target_bait.cytosine_info.Donor.CG.tsv
#bedtools intersect -a target_bait.cytosine_info.Donor.CG.bed -b /storage2/User/siyoo/module/demux2ready/src/MANE.GRCh38.v1.3.summary.chr.sorted.bed -wa -loj > target_bait.cytosine_info.Donor.CG.bed.anno
#awk 'BEGIN {FS="\t"; OFS="\t"; split("POU5F1 SOX2 KLF4 c-MYC TP53 BCL2L1", a, " "); for (i in a) b[a[i]]} $11 in b {print}' target_bait.cytosine_info.Donor.CG.bed.anno

#!/bin/bash

# Define the list of samples
samples=("Donor" "P05" "P10" "P15" "P20" "P25" "P30")

# Loop over each sample
for sample in "${samples[@]}"; do
    # Define the input file path based on the sample name
    input_file="../bismark/${sample}_split_cx_by_context/${sample}_R1_bismark_bt2_pe.deduplicated.bedCov.CG.gz"

    # Define the intermediate and final output file paths
    intermediate_output="target_bait.cytosine_info.${sample}.CG.bed"
    final_output="target_bait.cytosine_info.${sample}.CG.bed.anno.xls"

    # Run the bedtools map command
    bedtools map -a target_bait.bed -b "$input_file" -c 4,5,6 -o count,sum,sum | \
    awk 'BEGIN {FS="\t"; OFS="\t"} {depth = $5 + $6; if (depth > 0) {ratio = $5 / depth * 100} else {ratio = "NA"}; print $1, $2, $3, $4, $5, $6, ratio}' | \
    sed '1i#chr\tstart\tend\tcytosine_count\tmethyl_sum\tunmethyl_sum\tmethyl_level' > "$intermediate_output"

    # Run the bedtools intersect command
    bedtools intersect -a "$intermediate_output" -b /storage2/User/siyoo/module/demux2ready/src/MANE.GRCh38.v1.3.summary.chr.sorted.bed -wa -loj  | \
    sed '1i#chr\tstart\tend\tcytosine_count\tmethyl_sum\tunmethyl_sum\tmethyl_level\tgene_chr\tgene_start\tgene_end\tgene_name\tgene_score\tgene_strand' > "$final_output"

    # Filter the final output for specific gene names
    awk 'BEGIN {FS="\t"; OFS="\t"; split("POU5F1 SOX2 KLF4 c-MYC TP53 BCL2L1", a, " "); for (i in a) b[a[i]]} $11 in b {print}' "$final_output" | \
    sed '1i#chr\tstart\tend\tcytosine_count\tmethyl_sum\tunmethyl_sum\tmethyl_level\tgene_chr\tgene_start\tgene_end\tgene_name\tgene_score\tgene_strand' > "${final_output}.filtered.xls"

    echo "Processed $sample"
done
