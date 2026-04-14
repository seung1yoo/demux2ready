CONF="config.sh"
LOG="logs"
DR_dir="demux2ready"
RF_dir="readyfastq"

if [[ ! -f "$CONF" ]]; then
    echo ""
    echo "[INFO] config.sh not found. Running make_config.sh..."
    /bin/bash /home/rawfastq/demux2ready/workflow/make_config.sh
    echo ""
    echo "[INFO] Please edit $CONF (fill SeqType / Mode / Name / Affiliation / Date)."
    echo "[INFO] After editing, run /home/rawfastq/demux2ready/workflow/RawFastq.sh again."
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

case "${SeqType,,}" in
    paired* ) end=(R1 R2) ;;
    single* ) end=(R1) ;;
esac
ends="${end[*]}"

if [[ ! -d "${DR_dir}/${RF_dir}" ]]; then
    mkdir -p "${DR_dir}/${RF_dir}"
fi

if [[ ! -d "${LOG}" ]]; then
    mkdir -p $LOG
fi
echo "snakemake --executor cluster-generic --jobs 96 --cluster-generic-submit-cmd \"sbatch -J {rule} -c {threads} --open-mode=append -o $LOG/{rule}.out -e $LOG/{rule}.err\"" > Run.sh
chmod 755 Run.sh

##### Make Snakefile #####

cat > Snakefile <<EOF
import os, glob, re

PR = "${PR}"
raw_path = "${path}"
Oname = "${Oname}".split()
Rname = "${Rname}".split()
ends = "${ends}".split()
Mode = "${Mode}"
Name = "${Name}"
Affi = "${Affiliation}"
Date = "${Date}"
DR_dir = "${DR_dir}"
RF_dir = "${RF_dir}"
pairs = list(zip(Oname, Rname))

def find_sources(oname: str, end: str):
    paths = raw_path.split()
    idx = ("1", "R1") if end == "R1" else ("2", "R2")
    hits = []
    for rp in paths:
        for iidx in idx:
            patterns = [os.path.join(rp, "**", f"*{oname}*_{iidx}.fastq.gz"), os.path.join(rp, "**", f"*{oname}*_{iidx}.fq.gz"),]
            for p in patterns:
                hits.extend(glob.glob(p, recursive = True))
    unique_sorted = sorted(set(hits))
    if not unique_sorted:
        raise FileNotFoundError(f"No source FASTQs found for '{oname}'")
    return unique_sorted

r2o = {r: o for (o,r) in pairs}
def is_paired(r: str) -> bool:
    on = r2o[r]
    try:
        hits = find_sources(on, "R2")
        return len(hits) > 0
    except FileNotFoundError:
        return False

SOURCES = {}
for on, rn in pairs:
    SOURCES.setdefault(rn, {})
    for e in ends:
        SOURCES[rn][e] = find_sources(on, e)
REPORT = f"{PR}-RawFastq-QC-Report-v{Date}"
xlsx = f"{PR}-RawFastq-Download-Info-v{Date}.xlsx"

ruleorder: fastq_screen_paired > fastq_screen_single

rule all:
    input:
        [f"{DR_dir}/{RF_dir}/{r}_R1.fastq.gz" for r in Rname] + [f"{DR_dir}/{RF_dir}/{r}_R2.fastq.gz" for r in Rname if is_paired(r)],
#        [f"{DR_dir}/{RF_dir}/{r}_{e}.fastq.gz" for r in Rname for e in ends],
        f"{DR_dir}/fastp.summary.tsv",
        f"{DR_dir}/md5sum.summary.tsv",
        [f"{DR_dir}/fastq_screen/{r}_R1_screen.html" for r in Rname] + [f"{DR_dir}/fastq_screen/{r}_R2_screen.html" for r in Rname if is_paired(r)],
#        [f"{DR_dir}/fastq_screen/{r}_{e}_screen.html" for r in Rname for e in ends],
        f"{REPORT}.zip",
        f"{xlsx}"

rule copy_fastq:
    input:
        lambda wc: SOURCES[wc.r][wc.e]
    output:
        os.path.join(DR_dir, RF_dir, "{r}_{e}.fastq.gz")
    params:
        outdir = lambda wc: f"{DR_dir}/{RF_dir}"
    threads: 8
    run:
        ins = list(input)
        if len(ins) == 1:
            shell("cp -L {ins[0]} {output}")
        else:
            shell("cat {ins} > {output}")

rule fastp:
    input:
        lambda wc: [f"{DR_dir}/{RF_dir}/{wc.r}_R1.fastq.gz"] + ([f"{DR_dir}/{RF_dir}/{wc.r}_R2.fastq.gz"] if "R2" in ends else [])
    output:
        JSON = f"{DR_dir}/{RF_dir}" + "/{r}.fastp.json",
        HTML = f"{DR_dir}/{RF_dir}" + "/{r}.fastp.html"
    threads: 16
    run:
        ins = list(input)
        if len(ins) == 2:
            shell("fastp --thread {threads} -i {r1} --in2 {r2} --json {json} --html {html}".format(threads=threads, r1=ins[0], r2=ins[1], json=output.JSON, html=output.HTML))
        else:
            shell("fastp --thread {threads} -i {r1} --json {json} --html {html}".format(threads=threads, r1=ins[0], json=output.JSON, html=output.HTML))

rule fastp_parse:
    input:
        lambda wc: [f"{DR_dir}/{RF_dir}/{r}.fastp.json" for r in Rname]
    output:
        f"{DR_dir}/fastp.summary.tsv"
    params:
        mode = ("SE" if "R2" not in ends else "PE"),
        json = f"{DR_dir}/{RF_dir}",
        samples = " ".join(Rname),
        parse_sh = "/home/rawfastq/demux2ready/workflow/snakemake_fastp_parse.sh"
    threads: 1
    shell:
        "/bin/bash {params.parse_sh} {params.mode} {output} {params.json} {params.samples}"

rule md5sum:
    input:
        lambda wc: [f"{DR_dir}/{RF_dir}/{r}_{e}.fastq.gz" for r in Rname for e in ends]
    output:
        f"{DR_dir}/md5sum.summary.tsv"
    threads: 10
    shell:
        r"""
        rm -f {output}
        printf "%s\n" {input} | xargs -n 1 -P {threads} md5sum > {output}
        """

#rule fastq_screen:
#    input:
#        lambda wc: [f"{DR_dir}/{RF_dir}/{wc.r}_R1.fastq.gz"] + ([f"{DR_dir}/{RF_dir}/{wc.r}_R2.fastq.gz"] if "R2" in ends else [])
#    output:
#        html = f"{DR_dir}/fastq_screen" + "/{r}_{e}_screen.html"
#    threads: 16
#    run:
#        outdir = os.path.join(DR_dir, "fastq_screen")
#        ins = list(input)
#        if len(ins) == 2:
#            shell("fastq_screen --threads {threads} {r1} {r2} --outdir {outdir} --top 100000".format(threads=threads, r1=ins[0], r2=ins[1], outdir=outdir))
#        else:
#            shell("fastq_screen --threads {threads} {r1} --outdir {outdir} --top 100000".format(threads=threads, r1=ins[0], outdir=outdir))

rule fastq_screen_paired:
    input:
        r1 = f"{DR_dir}/{RF_dir}/{{r}}_R1.fastq.gz",
        r2 = f"{DR_dir}/{RF_dir}/{{r}}_R2.fastq.gz"
    output:
        r1_html = f"{DR_dir}/fastq_screen/{{r}}_R1_screen.html",
        r2_html = f"{DR_dir}/fastq_screen/{{r}}_R2_screen.html"
    threads: 16
    shell:
        "fastq_screen --threads {threads} {input.r1} {input.r2} --outdir {DR_dir}/fastq_screen --top 100000"

rule fastq_screen_single:
    input:
        r1 = f"{DR_dir}/{RF_dir}/{{r}}_R1.fastq.gz"
    output:
        r1_html = f"{DR_dir}/fastq_screen/{{r}}_R1_screen.html"
    threads: 16
    shell:
        "fastq_screen --threads {threads} {input.r1} --outdir {DR_dir}/fastq_screen --top 100000"

rule report:
    input:
        FASTP = f"{DR_dir}/fastp.summary.tsv",
        FASTQ_SCREEN = [f"{DR_dir}/fastq_screen/{r}_{e}_screen.html" for r in Rname for e in ends]
    output:
        ZIP = f"{REPORT}.zip",
        XLSX = f"{xlsx}"
    params:
        final_sh = "/home/rawfastq/demux2ready/workflow/snakemake_final.sh"
    threads: 8
    shell:
        "/bin/bash {params.final_sh}"
EOF

read -r -a RARR <<< "${Rname}"
echo -e "\n[OK] Snakefile generation completed!"
echo "Project: ${PR}"
echo "Samples: ${#RARR[@]} (${SeqType})"
echo -e "/bin/bash ./Run.sh\n"
