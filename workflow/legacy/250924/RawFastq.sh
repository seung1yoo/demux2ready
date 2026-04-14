CONF="config.sh"
LOG="rawfastq.log"
RF_dir="demux2ready/readyfastq"

if [[ ! -f "$CONF" ]]; then
    echo ""
    echo "[INFO] config.sh not found. Running 00.Config.sh..."
    /bin/bash /home/rawfastq/demux2ready/00.Config.sh
    echo ""
    echo "[INFO] Please edit $CONF (fill SeqType / Mode / Name / Affiliation / Date)."
    echo "[INFO] After editing, run /home/rawfastq/demux2ready/RawFastq.sh again."
    echo ""
    exit 0
fi

source "$CONF"

case "${SeqType,,}" in
    "paired-end"|"single-end") : ;;
    *)
    echo -e "\n[ERROR] SeqType must be either 'paired-end' or 'single-end'. Current: '${SeqType}'\n"
    exit 1
    ;;
esac
case "${Mode}" in
    "ngs-only"|"mb-consortium"|"dtc-microbiome") : ;;
    *)
    echo -e "\n[ERROR] Mode must be one of 'ngs-only', 'mb-consortium', or 'dtc-microbiome'. Current: '${Mode}'\n"
    exit 1
    ;;
esac
if [[ "${Mode}" != "mb-consortium" ]]; then
    if [[ -z "${Name:-}" ]]; then
        echo -e "\n[ERROR] Name is empty!!!\n"
        exit 1
    fi
    if [[ -z "${Affiliation:-}" ]]; then
        echo -e "\n[ERROR] Affiliation is empty!!!\n"
    exit 1
    fi
fi

read -r -a RARR <<< "${Rname}"
case "${SeqType,,}" in
    paired* ) end=(R1 R2) ;;
    single* ) end=(R1) ;;
esac

if [[ ! -d "$RF_dir" ]]; then
    mkdir -p "$RF_dir"
fi

check_missing(){
    local type="$1"
    missing=()
    for s in "${RARR[@]}"; do
        for e in "${end[@]}"; do
            f="${RF_dir}/${s}_${e}.${type}"
            if [[ ! -s "$f" ]]; then
                missing+=("$f")
            fi
        done
    done
}

expected=$((${#RARR[@]}*${#end[@]}))
check_missing fastq.gz
if [[ ! -f "cp_done" ]]; then
    if (( ${#missing[@]} == expected )); then
        echo -e "\n[COPY] Starting Copy... (Samples: ${#RARR[@]}, type: ${SeqType})\n"
    else
        echo -e "\n[COPY RETRY] Starting Copy... (Missing Files: ${#missing[@]}, Samples: ${#RARR[@]}, type: ${SeqType})\n"
    fi
    /bin/bash /home/rawfastq/demux2ready/01.Copy.sh
    if [[ ! -f "qc_done" ]];then
        echo -e "\n[QC] Starting QC... (Samples: ${#RARR[@]}, type: ${SeqType})\n"
        /bin/bash /home/rawfastq/demux2ready/02.QC.sh
        echo -e "\n[INFO] Check the QC results.\n[INFO] For REPORT generation, run /home/rawfastq/demux2ready/RawFastq.sh\n"
        exit 1
    fi
elif [[ ! -f "qc_done" ]]; then
    echo -e "\n[QC] Starting QC... (Samples: ${#RARR[@]}, type: ${SeqType})\n"
    /bin/bash /home/rawfastq/demux2ready/02.QC.sh
    echo -e "\n[INFO] Check the QC results.\n[INFO] For REPORT generation, run /home/rawfastq/demux2ready/RawFastq.sh\n"
    exit 1
fi

if [[ -f "qc_done" ]]; then
    echo -e "\n[REPORT] Starting report generation... (Samples: ${#RARR[@]}, type: ${SeqType})\n"
    /bin/bash /home/rawfastq/demux2ready/03.Report.sh
fi
