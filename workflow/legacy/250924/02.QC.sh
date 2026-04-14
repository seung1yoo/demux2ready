source config.sh

DR_dir="demux2ready"
RF_dir="readyfastq"
LOG="rawfastq.log"

read -r -a RARR <<< "${Rname}"

: > fastp.sh
if [[ "${SeqType,,}" == "single-end" ]]; then
    for i in "${!RARR[@]}"; do
        echo "if fastp --thread 16 -i ${DR_dir}/${RF_dir}/${RARR[$i]}_R1.fastq.gz --json ${DR_dir}/${RF_dir}/${RARR[$i]}.fastp.json --html ${DR_dir}/${RF_dir}/${RARR[$i]}.fastp.html; then echo FASTP_OK:${RARR[$i]}; else echo FASTP_FAILED:${RARR[$i]}; fi"
    done >> fastp.sh
else
    for i in "${!RARR[@]}"; do
        echo "if fastp --thread 16 -i ${DR_dir}/${RF_dir}/${RARR[$i]}_R1.fastq.gz --in2 ${DR_dir}/${RF_dir}/${RARR[$i]}_R2.fastq.gz --json ${DR_dir}/${RF_dir}/${RARR[$i]}.fastp.json --html ${DR_dir}/${RF_dir}/${RARR[$i]}.fastp.html; then echo FASTP_OK:${RARR[$i]}; else echo FASTP_FAILED:${RARR[$i]}; fi"
     done >> fastp.sh
fi

echo -e "\n==== $(date '+%F %T') | Fastp Start. Project: ${PR}, Samples: ${#RARR[@]} ${SeqType}" >> "$LOG"
echo "[QC-fastp] Starting fastp... (Samples: ${#RARR[@]}, type: ${SeqType})"
xargs -P 10 -L 1 -I CMD bash -c 'set -Eeuo pipefail; CMD' < fastp.sh | awk '{ print strftime("[%F %T]"), $0; fflush(); }' | tee -a "$LOG"

echo -e "\n[QC-fastp OK] fastp run finished in ${DR_dir}/${RF_dir}\n"
rm fastp.sh

echo -e "Samples\tbefore,ReadsCount\tbefore,BasesCount\tbefore,BasesCountGb\tbefore,GCRate\tbefore,Q30BaseRate\tbefore,MeanLen\tafter,ReadsCount\tafter,BasesCount\tafter,BasesCountGb\tafter,GCRate\tafter,Q30BaseRate\tafter,MeanLen\tDupRate\tInsertSize\tGoodReadRate" > "${DR_dir}/fastp.summary.tsv"
if [[ "${SeqType,,}" == "single-end" ]]; then
    for i in "${!RARR[@]}"; do
        j="${DR_dir}/${RF_dir}/${RARR[$i]}.fastp.json"
        jq -r --arg s "${RARR[$i]}" '[$s, .summary.before_filtering.total_reads, .summary.before_filtering.total_bases, ((.summary.before_filtering.total_bases/1e9*1000 | floor)/1000), ((.summary.before_filtering.gc_content*100*1000 | floor)/1000), ((.summary.before_filtering.q30_rate*100*1000 | floor)/1000), .summary.before_filtering.read1_mean_length, .summary.after_filtering.total_reads, .summary.after_filtering.total_bases, ((.summary.after_filtering.total_bases/1e9*1000 | floor)/1000), ((.summary.after_filtering.gc_content*100*1000 | floor)/1000), ((.summary.after_filtering.q30_rate*100*1000 | floor)/1000), .summary.after_filtering.read1_mean_length, ((.duplication.rate*100*1000 | floor)/1000), .insert_size.peak, ((.filtering_result.passed_filter_reads/.summary.before_filtering.total_reads*100*1000 | floor)/1000)] | @tsv' "$j"
    done >> "${DR_dir}/fastp.summary.tsv"
else
    for i in "${!RARR[@]}"; do
        j="${DR_dir}/${RF_dir}/${RARR[$i]}.fastp.json"
        bf=$(jq -r '[.summary.before_filtering.total_reads, .summary.before_filtering.total_bases, ((.summary.before_filtering.total_bases/1e9*1000 | floor)/1000), ((.summary.before_filtering.gc_content*100*1000 | floor)/1000), ((.summary.before_filtering.q30_rate*100*1000 | floor)/1000)] | @tsv' "$j")
        bfml=$(jq -r '[.summary.before_filtering.read1_mean_length, .summary.before_filtering.read2_mean_length] | @csv' "$j")
        af=$(jq -r '[.summary.after_filtering.total_reads, .summary.after_filtering.total_bases, ((.summary.after_filtering.total_bases/1e9*1000 | floor)/1000), ((.summary.after_filtering.gc_content*100*1000 | floor)/1000), ((.summary.after_filtering.q30_rate*100*1000 | floor)/1000)] | @tsv' "$j")
        afml=$(jq -r '[.summary.after_filtering.read1_mean_length, .summary.after_filtering.read2_mean_length] | @csv' "$j")
        etc=$(jq -r '[((.duplication.rate*100*1000 | floor)/1000), .insert_size.peak, ((.filtering_result.passed_filter_reads/.summary.before_filtering.total_reads*100*1000 | floor)/1000)] | @tsv' "$j")
        echo -e "${RARR[$i]}\t${bf}\t${bfml}\t${af}\t${afml}\t${etc}"
    done >> "${DR_dir}/fastp.summary.tsv"
fi

echo -e "\n==== $(date '+%F %T') | md5sum Start. Project: ${PR}, Samples: ${#RARR[@]} ${SeqType}" >> "$LOG"
echo "[QC-md5sum] Starting md5sum... (Samples: ${#RARR[@]}, type: ${SeqType})"
if [[ -f ${DR_dir}/md5sum.summary.tsv ]]; then
    rm ${DR_dir}/md5sum.summary.tsv
fi
find "${DR_dir}/${RF_dir}" -type f -name "*.fastq.gz" | xargs -n 1 -P 10 md5sum >> ${DR_dir}/md5sum.summary.tsv
echo -e "\n[QC-md5sum OK] MD5 checksum finished in ${DR_dir}/md5sum.summary.tsv\n"

: > fastq_screen.sh
if [[ "${SeqType,,}" == "single-end" ]]; then
    for i in "${!RARR[@]}"; do
        echo "if fastq_screen --thread 16 ${DR_dir}/${RF_dir}/${RARR[$i]}_R1.fastq.gz --outdir ${DR_dir}/fastq_screen --top 100000; then echo FQscreen_OK:${RARR[$i]}; else echo FQscreen_FAILED:${RARR[$i]}; fi"
    done >> fastq_screen.sh
else
    for i in "${!RARR[@]}"; do
        echo "if fastq_screen --thread 16 ${DR_dir}/${RF_dir}/${RARR[$i]}_R1.fastq.gz ${DR_dir}/${RF_dir}/${RARR[$i]}_R2.fastq.gz --outdir ${DR_dir}/fastq_screen --top 100000; then echo FQscreen_OK:${RARR[$i]}; else echo FQscreen_FAILED:${RARR[$i]}; fi"
    done >> fastq_screen.sh
fi

echo -e "\n==== $(date '+%F %T') | fastq_screen Start. Project: ${PR}, Samples: ${#RARR[@]} ${SeqType}" >> "$LOG"
echo "[QC-fastq_screen] Starting fastq_screen... (Samples: ${#RARR[@]}, type: ${SeqType})"
xargs -P 10 -L 1 -I CMD bash -c 'set -Eeuo pipefail; CMD' < fastq_screen.sh | awk '{ print strftime("[%F %T]"), $0; fflush(); }' | tee -a "$LOG"
touch qc_done
echo -e "\n[QC-fastq_screen OK] fastq_screen run finished in ${DR_dir}/fastq_screen\n"
rm fastq_screen.sh
