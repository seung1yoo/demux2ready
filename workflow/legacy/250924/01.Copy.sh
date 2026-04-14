source config.sh

RF_dir="demux2ready/readyfastq"
LOG="rawfastq.log"

read -r -a OARR <<< "${Oname}"
read -r -a RARR <<< "${Rname}"

: > copy.sh

if [[ "${SeqType,,}" == "single-end" ]]; then
    for i in "${!OARR[@]}"; do
        rd=$(find "$path" -type f -name "${OARR[$i]}_1.fastq.gz" -o -name "${OARR[$i]}_1.fq.gz" | sort)
        if [[ $(echo "$rd"|wc -w) -eq 1 ]]; then
            echo "if cp -vn $rd ${RF_dir}/${RARR[$i]}_R1.fastq.gz; then echo COPIED:${RARR[$i]}; else echo FAILED:${RARR[$i]}; fi"
        else
            echo "if cat $rd > ${RF_dir}/${RARR[$i]}_R1.fastq.gz; then echo COPIED:${RARR[$i]}; else echo FAILED:${RARR[$i]}; fi"
        fi
    done >> copy.sh
else
    for i in "${!OARR[@]}"; do
        rd1=$(find "$path" -type f -name "${OARR[$i]}_1.fastq.gz" -o -name "${OARR[$i]}_1.fq.gz" | sort)
        rd2=$(find "$path" -type f -name "${OARR[$i]}_2.fastq.gz" -o -name "${OARR[$i]}_2.fq.gz" | sort)
        if [[ $(echo "$rd1"|wc -w) -eq 1 ]]; then
            echo "if cp -vn $rd1 ${RF_dir}/${RARR[$i]}_R1.fastq.gz; then echo COPIED:${RARR[$i]}_R1; else echo FAILED:${RARR[$i]}_R1; fi"
        else
            echo "if cat $(echo $rd1 | tr '\n' ' ') > ${RF_dir}/${RARR[$i]}_R1.fastq.gz; then echo COPIED:${RARR[$i]}_R1; else echo FAILED:${RARR[$i]}_R1; fi"
        fi
        if [[ $(echo "$rd2"|wc -w) -eq 1 ]]; then
            echo "if cp -vn $rd2 ${RF_dir}/${RARR[$i]}_R2.fastq.gz; then echo COPIED:${RARR[$i]}_R2; else echo FAILED:${RARR[$i]}_R2; fi"
        else
            echo "if cat $(echo $rd2 | tr '\n' ' ') > ${RF_dir}/${RARR[$i]}_R2.fastq.gz; then echo COPIED:${RARR[$i]}_R2; else echo FAILED:${RARR[$i]}_R2; fi"
        fi
    done >> copy.sh
fi

echo "==== $(date '+%F %T') | Copy Start. Project: ${PR}, Samples: ${#OARR[@]} ${SeqType}" > "$LOG"
xargs -P 10 -L 1 -I CMD bash -c 'set -Eeuo pipefail; CMD' < copy.sh | awk '{ print strftime("[%F %T]"), $0; fflush(); }' | tee -a "$LOG"

case "${SeqType,,}" in
    paired* ) end=(R1 R2) ;;
    single* ) end=(R1) ;;
esac

expected=$((${#RARR[@]}*${#end[@]}))
present=()
missing=()
for s in "${RARR[@]}"; do
    for e in "${end[@]}"; do
        if [[ -s "${RF_dir}/${s}_${e}.fastq.gz" ]]; then
            ((present++))
        else
            missing+=("${s}_${e}")
        fi
    done
done
if (( present == expected )); then
    touch cp_done
    rm copy.sh
    echo -e "\n[COPY OK] All files copy finished in ${RF_dir} (Samples: ${#RARR[@]}, type: ${SeqType})\n"
else
    echo -e "\n[ERROR] Missing files detected (n: ${#missing[@]})"
    printf '  - %s\n' "${missing[@]}"
    echo ""
    exit 1
fi
