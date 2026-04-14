source config.sh

DR_dir="demux2ready"
RF_dir="readyfastq"
LOG="rawfastq.log"

OUT="${PR}-RawFastq-QC-Report-v${Date}"
if [[ -d "$OUT" ]]; then
    rm -r "$OUT"
fi
echo "==== $(date '+%F %T') | Report generation Start. Project: ${PR}, Samples: ${#OARR[@]} ${SeqType}" >> "$LOG"
if [[ "$Mode" == "mb-consortium" ]]; then
    python /home/rawfastq/atgcu-report/atgcu-report.py ${Mode} --outdir ${OUT} --project-id ${PR} --fastp-dir ${DR_dir}/${RF_dir} --fastq-screen-dir ${DR_dir}/fastq_screen
else
    python /home/rawfastq/atgcu-report/atgcu-report.py ${Mode} --outdir ${OUT} --project-id ${PR} --customer-name ${Name} --customer-affiliation ${Affiliation} --fastp-dir ${DR_dir}/${RF_dir} --fastq-screen-dir ${DR_dir}/fastq_screen
fi
zip -r ${OUT}.zip ${OUT}

pr=$(echo "${PR,,}")
id=$(id -un)
dt="20${Date}"
bucket="${pr}-${id}-${dt}"
echo "file_path,bucket_name" > atgcu-util.s3-upload.input.csv
for i in "${DR_dir}/${RF_dir}"/*.gz; do
    printf '%s,%s\n' "$i" "$bucket"
done >> atgcu-util.s3-upload.input.csv
echo "${OUT}.zip,${bucket}" >> atgcu-util.s3-upload.input.csv
python /home/rawfastq/atgcu-util/s3-upload/s3-upload.py --infn atgcu-util.s3-upload.input.csv

echo -e "\n[INFO] To refresh the presigned-URL expiration, please rerun:\npython /home/rawfastq/atgcu-util/s3-upload/s3-upload.py --infn atgcu-util.s3-upload.input.csv"
echo -e "\n[INFO] To generate the Download-Info.xlsx, please run:\npython /home/rawfastq/demux2ready/lib/create_excel_report.py --output ${PR}-RawFastq-Download-Info-v$(date +%y%m%d).xlsx"
