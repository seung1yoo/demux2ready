source config.sh

DR_dir="demux2ready"
RF_dir="readyfastq"
REPORT="${PR}-RawFastq-QC-Report-v${Date}"
xlsx="${PR}-RawFastq-Download-Info-v${Date}.xlsx"

if [[ -d "${REPORT}" ]]; then
    rm -r "${REPORT}"
fi

if [[ ${Mode} == "mb-consortium" ]]; then
    python /home/rawfastq/atgcu-report/atgcu-report.py ${Mode} --outdir "${REPORT}" --project-id ${PR} --fastp-dir ${DR_dir}/${RF_dir} --fastq-screen-dir ${DR_dir}/fastq_screen
else
    python /home/rawfastq/atgcu-report/atgcu-report.py ${Mode} --outdir "${REPORT}" --project-id ${PR} --customer-name ${Name} --customer-affiliation ${Affiliation} --fastp-dir ${DR_dir}/${RF_dir} --fastq-screen-dir ${DR_dir}/fastq_screen
fi

zip -r ${REPORT}.zip ${REPORT}

pr=$(echo "${PR,,}")
id=$(id -un)
dt="20${Date}"
bucket="${pr}-${id}-${dt}"
csv="atgcu-util.s3-upload.input.csv"
echo "file_path,bucket_name" > ${csv}
for i in "${DR_dir}/${RF_dir}"/*.gz; do
    printf '%s,%s\n' "$i" "$bucket" >> ${csv}
done
echo "${REPORT}.zip,${bucket}" >> ${csv}
python /home/rawfastq/atgcu-util/s3-upload/s3-upload.py --infn ${csv}
python /home/rawfastq/atgcu-util/s3-upload/s3-upload.py --infn ${csv}
python /home/rawfastq/demux2ready/lib/create_excel_report.py --output ${xlsx}
