MODE="${1:?missing MODE (SE|PE)}"
OUT="${2:?missing OUT TSV path}"
JSON="${3:?missing JSON dir}"
shift 3
SAMPLES=("$@")

printf "Samples\tReadsCount\tBasesCount\tBasesCountGb\tGCRate\tQ20BaseRate\tQ30BaseRate\tDupRate\tInsertSize\tR1MeanLen\tR2MeanLen\tGoodReadRate\n" > "$OUT"

if [[ "$MODE" == "SE" ]]; then
    for s in "${SAMPLES[@]}"; do
        j="${JSON}/${s}.fastp.json"
        jq -r --arg s "$s" '[$s, .summary.before_filtering.total_reads, .summary.before_filtering.total_bases, ((.summary.before_filtering.total_bases/1e9*1000 | floor)/1000), ((.summary.before_filtering.gc_content*100*1000 | floor)/1000), ((.summary.before_filtering.q20_rate*100*1000 | floor)/1000), ((.summary.before_filtering.q30_rate*100*1000 | floor)/1000), ((.duplication.rate*100*1000 | floor)/1000), .insert_size.peak, .summary.before_filtering.read1_mean_length, .summary.before_filtering.read1_mean_length, ((.filtering_result.passed_filter_reads/.summary.before_filtering.total_reads*100*1000 | floor)/1000)] | @tsv' "$j" >> "$OUT"
    done
else
    for s in "${SAMPLES[@]}"; do
        j="${JSON}/${s}.fastp.json"
        jq -r --arg s "$s" '[$s, .summary.before_filtering.total_reads, .summary.before_filtering.total_bases, ((.summary.before_filtering.total_bases/1e9*1000 | floor)/1000), ((.summary.before_filtering.gc_content*100*1000 | floor)/1000), ((.summary.before_filtering.q20_rate*100*1000 | floor)/1000), ((.summary.before_filtering.q30_rate*100*1000 | floor)/1000), ((.duplication.rate*100*1000 | floor)/1000), .insert_size.peak, .summary.before_filtering.read1_mean_length, .summary.before_filtering.read2_mean_length, ((.filtering_result.passed_filter_reads/.summary.before_filtering.total_reads*100*1000 | floor)/1000)] | @tsv' "$j" >> "$OUT"
    done
fi
