
#demux_fastq = demuxed fastq
#ready_fastq = raw fastq for pipeline
#ssheet_id = a ID in SampleSheet.csv
#target_id = a ID to converting
#DEMUXPATH = demux path
#IDMAP = ssheet_id:target_id
#READYPATH = ready path

#differentialabundance = a.k.a diffabdc

import subprocess
from pathlib import Path
import yaml
import json
import csv
import logging
logging.basicConfig(level=logging.DEBUG)

class Demux2Ready:
    def __init__(self, conf_fn, work_dir, find_fastq_mode):
        self.find_fastq_mode = find_fastq_mode
        Path(work_dir).mkdir(exist_ok=True)
        self.readypath = Path(work_dir) / "readyfastq"
        self.readypath.mkdir(exist_ok=True)

        # make a shell script for demux2ready
        self.demux2ready_sh = Path(work_dir) / "demux2ready.sh"
        # make a shell script for fastp statistics
        self.fastp_sh = Path(work_dir) / "fastp.sh"
        # sum-up fastp json files into A TSV file
        self.fastp_summary = Path(work_dir) / "fastp.summary.tsv"
        # make a shell script for fastq-screen
        self.fastqscreenpath = Path(work_dir) / "fastq_screen"
        self.fastqscreen_sh = Path(work_dir) / "fastq_screen.sh"
        # md5sum
        self.md5sum_summary = Path(work_dir) / "md5sum.summary.tsv"

        # make a config file for rnaseq-pipeline
        self.samplefile_path = Path(work_dir) / "sample_file"
        # make a config file for nf-core/rnaseq
        self.nfcore_rnaseq_samplesheet = Path(work_dir) / "nfcore_rnaseq_samplesheet.csv"
        # make a config file for nf-core/differentialabundance
        self.nfcore_diffabdc_samplesheet = Path(work_dir) / "nfcore_diffabdc_samplesheet.csv"
        self.nfcore_diffabdc_contrasts = Path(work_dir) / "nfcore_diffabdc_contrasts.csv"
        # make a config file for atgcu-rnaseq
        self.atgcu_rnaseq_sample = Path(work_dir) / "atgcu_rnaseq_sample.yml"
        self.atgcu_rnaseq_compare = Path(work_dir) / "atgcu_rnaseq_compare.yml"

        # make shell script for sequencing throughtput adjust
        self.baseadjustpath = Path(work_dir) / "baseadjustfastq"
        self.splitpath = self.baseadjustpath / "splitfastq"
        self.fastp_split_sh = Path(work_dir) / "fastp_split.sh"
        self.baseadjust_sh = Path(work_dir) / "baseadjust.sh"
        # make fastp shell script for read pre-processing
        self.cleanpath = Path(work_dir) / "cleanfastq"
        self.fastp_clean_sh = Path(work_dir) / "fastp_clean.sh"
        self.cutadapt_sh = Path(work_dir) / "cutadapt.sh"
        # make tophat shell script
        self.tophat_sh = Path(work_dir) / "tophat.sh"
        self.tophatpath = Path(work_dir) / "tophat"
        # make cuffquant shell script
        self.cuffquant_sh = Path(work_dir) / "cuffquant.sh"
        self.cuffquantpath = Path(work_dir) / "cuffquant"
        # make cufflinks shell script
        self.cufflinks_sh = Path(work_dir) / "cufflinks.sh"
        self.cufflinkspath = Path(work_dir) / "cufflinks"
        # make cuffnorm shell script
        self.cuffnorm_sh = Path(work_dir) / "cuffnorm.sh"
        self.cuffnormpath = Path(work_dir) / "cuffnorm"
        # make cuffdiff shell script
        self.cuffdiff_sh = Path(work_dir) / "cuffdiff.sh"
        self.cuffdiffpath = Path(work_dir) / "cuffdiff"
        # make bismark shell script
        self.bismarkpath = Path(work_dir) / "bismark"
        self.bismark_sh = Path(work_dir) / "bismark.sh"

        #
        self.demuxpath_s = self.grep_demuxpath(conf_fn)
        self.idmap_dic = self.grep_idmap(conf_fn)
        self.fastq_dic = self.init_fastq_dic()
        if args.mode in ['mksh_d2r','prepare']:
            self.find_fastq()
        self.comp_dic = self.grep_comp(conf_fn)
        self.meta_dic = self.grep_meta(conf_fn)
        self.baseadjust_dic = self.grep_baseadjust(conf_fn)

    def init_fastq_dic(self):
        fastq_dic = dict()
        fastq_dic.setdefault('demux', {})
        fastq_dic.setdefault('ready', {})
        fastq_dic.setdefault('split', {})
        fastq_dic.setdefault('baseadjust', {})
        fastq_dic.setdefault('clean', {})
        for ssheet_id, target_id in self.idmap_dic['for'].items():
            fastq_dic['demux'].setdefault(ssheet_id, {}).setdefault('r1_s', [])
            fastq_dic['demux'].setdefault(ssheet_id, {}).setdefault('r2_s', [])
            #
            _fastq_r1 = self.readypath / f"{target_id}_R1.fastq.gz"
            _fastq_r2 = self.readypath / f"{target_id}_R2.fastq.gz"
            fastq_dic['ready'].setdefault(target_id, {}).setdefault('r1', _fastq_r1)
            fastq_dic['ready'].setdefault(target_id, {}).setdefault('r2', _fastq_r2)
            #
            for i in range(10): # 10 is a default value
                _fastq_r1 = self.splitpath / f"{i+1:0>4}.{target_id}_R1.fastq.gz"
                _fastq_r2 = self.splitpath / f"{i+1:0>4}.{target_id}_R2.fastq.gz"
                fastq_dic['split'].setdefault(target_id, {}).setdefault(i+1, {}).setdefault('r1', _fastq_r1)
                fastq_dic['split'].setdefault(target_id, {}).setdefault(i+1, {}).setdefault('r2', _fastq_r2)
            #
            _fastq_r1 = self.baseadjustpath / f"{target_id}_R1.fastq.gz"
            _fastq_r2 = self.baseadjustpath / f"{target_id}_R2.fastq.gz"
            fastq_dic['baseadjust'].setdefault(target_id, {}).setdefault('r1', _fastq_r1)
            fastq_dic['baseadjust'].setdefault(target_id, {}).setdefault('r2', _fastq_r2)
            #
            _fastq_r1 = self.cleanpath / f"{target_id}_R1.fastq.gz"
            _fastq_r2 = self.cleanpath / f"{target_id}_R2.fastq.gz"
            fastq_dic['clean'].setdefault(target_id, {}).setdefault('r1', _fastq_r1)
            fastq_dic['clean'].setdefault(target_id, {}).setdefault('r2', _fastq_r2)
        return fastq_dic

    def find_fastq(self):
        for demuxpath in self.demuxpath_s:
            logging.debug(f"{demuxpath}")
            for ssheet_id, target_id in self.idmap_dic['for'].items():
                # for read
                if self.find_fastq_mode in ['outsourcing']:
                    r1_s = list(demuxpath.rglob(f'{ssheet_id}_R1.fastq.gz'))
                    r1_s.extend(list(demuxpath.rglob(f'{ssheet_id}_R1.fq.gz')))
                    r1_s.extend(list(demuxpath.rglob(f'{ssheet_id}_1.fastq.gz')))
                    r1_s.extend(list(demuxpath.rglob(f'{ssheet_id}_1.fq.gz')))
                elif self.find_fastq_mode in ['demux']:
                    r1_s = list(demuxpath.rglob(f'{ssheet_id}_S*_L*_R1_001.fastq.gz'))
                    r1_s.extend(list(demuxpath.rglob(f'{ssheet_id}_S*_R1_001.fastq.gz')))
                self.add_fastq_path_s('r1', r1_s, ssheet_id)

                # rev read
                if self.find_fastq_mode in ['outsourcing']:
                    r2_s = list(demuxpath.rglob(f'{ssheet_id}_R2.fastq.gz'))
                    r2_s.extend(list(demuxpath.rglob(f'{ssheet_id}_R2.fq.gz')))
                    r2_s.extend(list(demuxpath.rglob(f'{ssheet_id}_2.fastq.gz')))
                    r2_s.extend(list(demuxpath.rglob(f'{ssheet_id}_2.fq.gz')))
                elif self.find_fastq_mode in ['demux']:
                    r2_s = list(demuxpath.rglob(f'{ssheet_id}_S*_L*_R2_001.fastq.gz'))
                    r2_s.extend(list(demuxpath.rglob(f'{ssheet_id}_S*_R2_001.fastq.gz')))
                self.add_fastq_path_s('r2', r2_s, ssheet_id)

    def add_fastq_path_s(self, r_tag, r_s, ssheet_id):
        for r_path in r_s:
            logging.debug(f"{ssheet_id}, {r_tag}, {r_path}")
            if r_tag in ['r1']:
                if r_path not in self.fastq_dic['demux'][ssheet_id]['r1_s']:
                    self.fastq_dic['demux'][ssheet_id]['r1_s'].append(r_path)
            elif r_tag in ['r2']:
                if r_path not in self.fastq_dic['demux'][ssheet_id]['r2_s']:
                    self.fastq_dic['demux'][ssheet_id]['r2_s'].append(r_path)
            else:
                pass
        return 1

    def grep_demuxpath(self, conf_fn):
        demuxpath_s = list()
        for line in open(conf_fn):
            items = line.rstrip().split()
            if not items:
                continue
            if not items[0] in ['DEMUXPATH']:
                continue
            if items[1] not in demuxpath_s:
                demuxpath = Path(items[1])
                demuxpath_s.append(demuxpath)
        return demuxpath_s

    def grep_idmap(self, conf_fn):
        idmap_dic = dict()
        idmap_dic.setdefault('for', {})
        idmap_dic.setdefault('rev', {})
        for line in open(conf_fn):
            items = line.rstrip().split()
            if not items:
                continue
            if not items[0] in ['IDMAP']:
                continue
            idmap_dic['for'].setdefault(items[1], items[2])
            idmap_dic['rev'].setdefault(items[2], []).append(items[1])
        return idmap_dic

    def grep_comp(self, conf_fn):
        comp_dic = dict()

        for line in open(conf_fn):
            items = line.rstrip().split()
            if not items:
                continue
            if not items[0] in ['COMP']:
                continue

            comp_id = items[1]
            cont_group_id = items[2]
            cont_sample_s = items[3].split(',')
            case_group_id = items[4]
            case_sample_s = items[5].split(',')

            comp_dic.setdefault(comp_id, {}).setdefault('cont', {}).setdefault('group_id', cont_group_id)
            comp_dic.setdefault(comp_id, {}).setdefault('cont', {}).setdefault('sample_s', cont_sample_s)
            comp_dic.setdefault(comp_id, {}).setdefault('case', {}).setdefault('group_id', case_group_id)
            comp_dic.setdefault(comp_id, {}).setdefault('case', {}).setdefault('sample_s', case_sample_s)
        return comp_dic

    def grep_meta(self, conf_fn):
        meta_dic = dict()
        meta_dic.setdefault("cpu", "[CPU]")
        meta_dic.setdefault("gtf", "[GTF]")
        meta_dic.setdefault("mask_gtf", "[MASK_GTF]")
        meta_dic.setdefault("bowtie2_index", "[BOWTIE2_INDEX]")
        meta_dic.setdefault("genome_fasta", "[GENOME_FASTA]")
        meta_dic.setdefault("library_type", "[LIBRARY_TYPE]")
        meta_dic.setdefault("bismark_genome", "[BISMARK_GENOME]")
        meta_dic.setdefault("interval_list", "[INTERVAL_LIST]")
        for line in open(conf_fn):
            items = line.rstrip().split()
            if not items:
                continue
            if not items[0] in ['META']:
                continue
            key = items[1]
            if key in ["CPU"]:
                meta_dic["cpu"] = items[2]
            elif key in ["GTF"]:
                meta_dic["gtf"] = items[2]
            elif key in ["MASK_GTF"]:
                meta_dic["mask_gtf"] = items[2]
            elif key in ["BOWTIE2_INDEX"]:
                meta_dic["bowtie2_index"] = items[2]
            elif key in ["GENOME_FASTA"]:
                meta_dic["genome_fasta"] = items[2]
            elif key in ["LIBRARY_TYPE"]:
                meta_dic["library_type"] = items[2]
            elif key in ["BISMARK_GENOME"]:
                meta_dic["bismark_genome"] = items[2]
            elif key in ["INTERVAL_LIST"]:
                meta_dic["interval_list"] = items[2]
            elif key in ["TARGET_BED"]:
                meta_dic["target_bed"] = items[2]
        return meta_dic

    def grep_baseadjust(self, conf_fn):
        _dic = dict()
        for line in open(conf_fn):
            items = line.rstrip().split()
            if not items:
                continue
            if not items[0] in ['BASEADJUST']:
                continue
            _id = items[1]
            need_split_n = int(items[2])
            _dic.setdefault(_id, need_split_n)
        return _dic

    def make_demux2ready_sh(self):
        outfh = self.demux2ready_sh.open('w')
        for target_id, read_dic in self.fastq_dic['ready'].items():
            r1_s = list()
            r2_s = list()
            for ssheet_id in self.idmap_dic['rev'][target_id]:
                r1_s.extend(self.fastq_dic['demux'][ssheet_id]['r1_s'])
                r2_s.extend(self.fastq_dic['demux'][ssheet_id]['r2_s'])

            # for fastq read 1
            if len(r1_s) in [1]:
                _cmd = ['ln -s']
                _cmd.append(str(r1_s[0]))
                _cmd.append(str(read_dic['r1']))
                outfh.write(f"{' '.join(_cmd)}\n")
            else:
                _cmd = ['cat']
                _cmd.append(' '.join([str(x) for x in sorted(r1_s)]))
                _cmd.append('>')
                _cmd.append(str(read_dic['r1']))
                outfh.write(f"{' '.join(_cmd)}\n")
            # for fastq read 2
            if len(r2_s) in [1]:
                _cmd = ['ln -s']
                _cmd.append(str(r2_s[0]))
                _cmd.append(str(read_dic['r2']))
                outfh.write(f"{' '.join(_cmd)}\n")
            else:
                _cmd = ['cat']
                _cmd.append(' '.join([str(x) for x in sorted(r2_s)]))
                _cmd.append('>')
                _cmd.append(str(read_dic['r2']))
                outfh.write(f"{' '.join(_cmd)}\n")
        outfh.close()
        return 1

    def prepare_now(self):
        for target_id, read_dic in self.fastq_dic['ready'].items():
            r1_s = list()
            r2_s = list()
            for ssheet_id in self.idmap_dic['rev'][target_id]:
                r1_s.extend(self.fastq_dic['demux'][ssheet_id]['r1_s'])
                r2_s.extend(self.fastq_dic['demux'][ssheet_id]['r2_s'])

            # for fastq read 1
            if len(r1_s) in [1]:
                read_dic['r1'].symlink_to(r1_s[0])
            else:
                pass # to do using subprocess
            # for fastq read 2
            if len(r2_s) in [1]:
                read_dic['r2'].symlink_to(r2_s[0])
            else:
                pass # to do using subprocess
        return 1

    def make_samplefile(self):
        outfh = self.samplefile_path.open('w')
        outfh.write('REF\tH_sapiens_GRCh38_ENS109\n')
        outfh.write('MODE\tquant\n')
        sample_num = 0
        for target_id, read_dic in self.fastq_dic['ready'].items():
            sample_num += 1
            items = ['SAMPLE']
            items.append("IBC{0:0>4}".format(sample_num))
            items.append('truseq')
            items.append(target_id)
            items.append(target_id.split('_')[0])
            items.append(str(sample_num))
            items.append('1')
            items.append('1')
            items.append(str(read_dic['r1']))
            outfh.write("{0}\n".format('\t'.join(items)))
            items = ['SAMPLE']
            items.append("IBC{0:0>4}".format(sample_num))
            items.append('truseq')
            items.append(target_id)
            items.append(target_id.split('_')[0])
            items.append(str(sample_num))
            items.append('1')
            items.append('2')
            items.append(str(read_dic['r2']))
            outfh.write("{0}\n".format('\t'.join(items)))
        outfh.close()
        return 1

    def make_nfcore_rnaseq_samplesheet(self):
        outfh = self.nfcore_rnaseq_samplesheet.open("w")
        csvh = csv.writer(outfh)
        header = ["sample","fastq_1","fastq_2","strandedness"]
        csvh.writerow(header)
        for target_id, read_dic in self.fastq_dic['ready'].items():
            items = [target_id]
            items.append(str(read_dic['r1']))
            items.append(str(read_dic['r2']))
            items.append('auto')
            csvh.writerow(items)
        outfh.close()

    def make_nfcore_diffabdc_samplesheet(self):
        outfh = self.nfcore_diffabdc_samplesheet.open("w")
        csvh = csv.writer(outfh)
        header = ["sample","fastq_1","fastq_2","condition","replicate"]
        csvh.writerow(header)
        for target_id, read_dic in self.fastq_dic['ready'].items():
            items = [target_id]
            items.append(str(read_dic['r1']))
            items.append(str(read_dic['r2']))
            items.append('[control/case/...]')
            items.append('[1/2/3/...]')
            csvh.writerow(items)
        outfh.close()

    def make_nfcore_diffabdc_contrasts(self):
        outfh = self.nfcore_diffabdc_contrasts.open("w")
        csvh = csv.writer(outfh)
        header = ["id","variable","reference","target"]
        csvh.writerow(header)
        csvh.writerow(["DEG001","condition","control","case"])
        outfh.close()

    def make_atgcu_rnaseq_sample(self):
        outfh = self.atgcu_rnaseq_sample.open("w")
        _dic = dict()
        _dic.setdefault("samples", {})
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _dic["samples"].setdefault(target_id, {}).setdefault("fastq_1", str(read_dic['r1']))
            _dic["samples"].setdefault(target_id, {}).setdefault("fastq_2", str(read_dic['r2']))
            _dic["samples"].setdefault(target_id, {}).setdefault("cuff_librarytype", "fr-unstranded")
        yaml.dump(_dic, outfh)
        outfh.close()

    def make_atgcu_rnaseq_compare(self):
        outfh = self.atgcu_rnaseq_compare.open("w")
        _dic = dict()
        _dic.setdefault("compares", {})
        for comp_id, group_dic in self.comp_dic.items():
            _dic["compares"].setdefault(comp_id, {}).setdefault("control", {})
            _dic["compares"][comp_id]["control"].setdefault("group_name", group_dic["cont"]["group_id"])
            _dic["compares"][comp_id]["control"].setdefault("samples", group_dic["cont"]["sample_s"])
            _dic["compares"].setdefault(comp_id, {}).setdefault("treatment", {})
            _dic["compares"][comp_id]["treatment"].setdefault("group_name", group_dic["case"]["group_id"])
            _dic["compares"][comp_id]["treatment"].setdefault("samples", group_dic["case"]["sample_s"])
        yaml.dump(_dic, outfh)
        outfh.close()

    def make_fastp_sh(self):
        outfh = self.fastp_sh.open('w')
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['fastp']
            _cmd.append('--thread')
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append('-i')
            _cmd.append(str(read_dic['r1']))
            _cmd.append('--in2')
            _cmd.append(str(read_dic['r2']))
            _cmd.append('--json')
            _cmd.append(str(self.readypath / f"{target_id}.fastp.json"))
            _cmd.append('--html')
            _cmd.append(str(self.readypath / f"{target_id}.fastp.html"))
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_fastq_screen_sh(self):
        outfh = self.fastqscreen_sh.open('w')
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['fastq_screen']
            _cmd.append('--thread')
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append(str(read_dic['r1']))
            _cmd.append(str(read_dic['r2']))
            _cmd.append('--outdir')
            _cmd.append(str(self.fastqscreenpath))
            _cmd.append('--top 100000')
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_cutadapt_sh(self):
        outfh = self.cutadapt_sh.open('w')
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['cutadapt']
            _cmd.append('--cores')
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append('-u')
            _cmd.append('-119') # 151 PE to 32 PE
            _cmd.append('-U')
            _cmd.append('-119') # 151 PE to 32 PE
            _cmd.append('-o')
            _cmd.append(str(self.cleanpath / f"{target_id}_R1.fastq.gz"))
            _cmd.append('-p')
            _cmd.append(str(self.cleanpath / f"{target_id}_R2.fastq.gz"))
            _cmd.append(str(read_dic['r1']))
            _cmd.append(str(read_dic['r2']))
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_fastp_clean_sh(self):
        outfh = self.fastp_clean_sh.open('w')
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['fastp']
            _cmd.append('--thread')
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append('-i')
            _cmd.append(str(read_dic['r1']))
            _cmd.append('--out1')
            _cmd.append(str(self.cleanpath / f"{target_id}_R1.fastq.gz"))
            _cmd.append('--in2')
            _cmd.append(str(read_dic['r2']))
            _cmd.append('--out2')
            _cmd.append(str(self.cleanpath / f"{target_id}_R2.fastq.gz"))
            _cmd.append('--json')
            _cmd.append(str(self.cleanpath / f"{target_id}.fastp.json"))
            _cmd.append('--html')
            _cmd.append(str(self.cleanpath / f"{target_id}.fastp.html"))
            #_cmd.append('--trim_tail1 50') # optional read length adjust
            #_cmd.append('--trim_tail2 50') # optional read length adjust
            #_cmd.append('--disable_adapter_trimming') # optional read length adjust
            #_cmd.append('--disable_trim_poly_g') # optional read length adjust
            #_cmd.append('--disable_quality_filtering') # optional read length adjust
            #_cmd.append('--disable_length_filtering') # optional read length adjust
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_fastp_split_sh(self, split_n = 10): # 10 is default value. It is related to fastq_dic['split']
        outfh = self.fastp_split_sh.open('w')
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['fastp']
            _cmd.append('--thread')
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append('-i')
            _cmd.append(str(read_dic['r1']))
            _cmd.append('--out1')
            _cmd.append(str(self.splitpath / f"{target_id}_R1.fastq.gz"))
            _cmd.append('--in2')
            _cmd.append(str(read_dic['r2']))
            _cmd.append('--out2')
            _cmd.append(str(self.splitpath / f"{target_id}_R2.fastq.gz"))
            _cmd.append('--json')
            _cmd.append(str(self.splitpath / f"{target_id}.fastp.json"))
            _cmd.append('--html')
            _cmd.append(str(self.splitpath / f"{target_id}.fastp.html"))
            _cmd.append(f'--split {split_n}')
            _cmd.append('--disable_adapter_trimming')
            _cmd.append('--disable_trim_poly_g')
            _cmd.append('--disable_quality_filtering')
            _cmd.append('--disable_length_filtering')
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_baseadjust_sh(self):
        outfh = self.baseadjust_sh.open('w')
        for target_id, read_dic in self.fastq_dic['split'].items():
            need_split_n = self.baseadjust_dic[target_id]
            for r in ['r1', 'r2']:
                _cmd = ["cat"]
                for i in range(need_split_n):
                    _cmd.append(str(read_dic[i+1][r]))
                _cmd.append(">")
                _cmd.append(str(self.fastq_dic['baseadjust'][target_id][r]))
                outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def parse_fastp(self):
        fastp_dic = dict()
        for target_id, read_dic in self.fastq_dic['ready'].items():
            fastp_json = self.readypath / f"{target_id}.fastp.json"
            if fastp_json.exists():
                fastp_dic.setdefault(target_id, json.load(fastp_json.open("r")))
            fastp_json = self.cleanpath / f"{target_id}.fastp.json"
            if fastp_json.exists():
                fastp_dic.setdefault(target_id, json.load(fastp_json.open("r")))

        outfh = self.fastp_summary.open('w')
        headers = ["Samples"]
        headers.append("ReadsCount")
        headers.append("BasesCount")
        headers.append("BasesCountGb")
        headers.append("GCRate")
        headers.append("Q20BaseRate")
        headers.append("Q30BaseRate")
        headers.append("DupRate")
        headers.append("InsertSize")
        headers.append("R1MeanLen")
        headers.append("R2MeanLen")
        headers.append("GoodReadRate")
        outfh.write('{0}\n'.format("\t".join(headers)))
        for target_id, info_dic in fastp_dic.items():
            items = [target_id]
            total_reads = info_dic["summary"]["before_filtering"]["total_reads"]
            items.append(total_reads)
            items.append(info_dic["summary"]["before_filtering"]["total_bases"])
            items.append(round(int(info_dic["summary"]["before_filtering"]["total_bases"])*0.000000001, 3))
            items.append(round(float(info_dic["summary"]["before_filtering"]["gc_content"])*100, 3))
            items.append(round(float(info_dic["summary"]["before_filtering"]["q20_rate"])*100, 3))
            items.append(round(float(info_dic["summary"]["before_filtering"]["q30_rate"])*100, 3))
            items.append(round(float(info_dic["duplication"]["rate"])*100, 3))
            items.append(info_dic["insert_size"]["peak"])
            items.append(info_dic["summary"]["before_filtering"]["read1_mean_length"])
            items.append(info_dic["summary"]["before_filtering"]["read2_mean_length"])
            items.append(round(int(info_dic["filtering_result"]["passed_filter_reads"])/float(total_reads)*100, 3))
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()


    def make_tophat_sh(self):
        outfh = self.tophat_sh.open('w')
        for target_id, read_dic in self.fastq_dic['ready'].items():
            wkdir = self.tophatpath / target_id
            wkdir.mkdir(exist_ok=True)
            _cmd = ['tophat']
            _cmd.append("-o")
            _cmd.append(str(wkdir))
            _cmd.append("-p")
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append("--library-type")
            _cmd.append(self.meta_dic["library_type"]) # fr-unstranded, fr-firststrand, fr-secondstrand
            _cmd.append("--rg-id")
            _cmd.append(target_id)
            _cmd.append("--rg-sample")
            _cmd.append(target_id)
            _cmd.append("-G")
            _cmd.append(self.meta_dic["gtf"])
            _cmd.append(self.meta_dic["bowtie2_index"])
            _cmd.append(str(self.cleanpath / f"{target_id}_R1.fastq.gz"))
            _cmd.append(str(self.cleanpath / f"{target_id}_R2.fastq.gz"))
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_cuffquant_sh(self):
        outfh = self.cuffquant_sh.open('w')
        for target_id, read_dic in self.fastq_dic['ready'].items():
            wkdir = self.cuffquantpath / target_id
            wkdir.mkdir(exist_ok=True)
            _cmd = ['cuffquant']
            _cmd.append("-o")
            _cmd.append(str(wkdir))
            _cmd.append("-p")
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append("--library-type")
            _cmd.append(self.meta_dic["library_type"])
            _cmd.append("--multi-read-correct")
            _cmd.append("--frag-bias-correct")
            _cmd.append(self.meta_dic["genome_fasta"])
            _cmd.append("--mask-file")
            _cmd.append(self.meta_dic["mask_gtf"])
            _cmd.append(self.meta_dic["gtf"])
            _cmd.append(str(self.tophatpath / target_id / "accepted_hits.bam"))
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_cufflinks_sh(self):
        outfh = self.cufflinks_sh.open('w')
        for target_id, read_dic in self.fastq_dic['ready'].items():
            wkdir = self.cufflinkspath / target_id
            wkdir.mkdir(exist_ok=True)
            _cmd = ['cufflinks']
            _cmd.append("-o")
            _cmd.append(str(wkdir))
            _cmd.append("-p")
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append("--library-type")
            _cmd.append(self.meta_dic["library_type"])
            _cmd.append("--multi-read-correct")
            _cmd.append("--frag-bias-correct")
            _cmd.append(self.meta_dic["genome_fasta"])
            _cmd.append("--mask-file")
            _cmd.append(self.meta_dic["mask_gtf"])
            _cmd.append("-G")
            _cmd.append(self.meta_dic["gtf"])
            _cmd.append(str(self.tophatpath / target_id / "accepted_hits.bam"))
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_cuffnorm_sh(self):
        wkdir = self.cuffnormpath
        wkdir.mkdir(exist_ok=True)

        target_id_s = list()
        cxb_s = list()
        for target_id, read_dic in self.fastq_dic['ready'].items():
            target_id_s.append(target_id)
            cxb_s.append(str(self.cuffquantpath / target_id / "abundances.cxb"))
            #cxb_s.append(str(self.tophatpath / target_id / "accepted_hits.bam"))

        _cmd = ['cuffnorm']
        _cmd.append("-o")
        _cmd.append(str(wkdir))
        _cmd.append("-p")
        _cmd.append(self.meta_dic["cpu"])
        _cmd.append("--library-type")
        _cmd.append(self.meta_dic["library_type"])
        _cmd.append("--library-norm-method")
        _cmd.append("classic-fpkm")
        _cmd.append("--output-format")
        _cmd.append("simple-table")
        _cmd.append("--labels")
        _cmd.append(",".join(target_id_s))
        _cmd.append(self.meta_dic["gtf"])
        _cmd.extend(cxb_s)

        outfh = self.cuffnorm_sh.open('w')
        outfh.write("{0}\n".format(" ".join(_cmd)))
        outfh.close()



    def make_cuffdiff_sh(self):
        outfh = self.cuffdiff_sh.open("w")
        for comp_id, group_dic in self.comp_dic.items():
            wkdir = self.cuffdiffpath / comp_id
            wkdir.mkdir(exist_ok=True)
            _cmd = ['cuffdiff']
            _cmd.append("-o")
            _cmd.append(str(wkdir))
            _cmd.append("-p")
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append("--library-type")
            _cmd.append(self.meta_dic["library_type"])
            #_cmd.append("--mask-file")
            #_cmd.append(self.meta_dic["mask_gtf"])
            _cmd.append("--dispersion-method")
            _cmd.append("pooled")
            _cmd.append("--library-norm-method")
            _cmd.append("geometric")
            _cmd.append("-L")
            cont_group_id = group_dic["cont"]["group_id"]
            case_group_id = group_dic["case"]["group_id"]
            _cmd.append(f"{cont_group_id},{case_group_id}")
            _cmd.append(self.meta_dic["gtf"])

            cont_bam = ','.join([str(self.cuffquantpath / x / "abundances.cxb") for x in group_dic["cont"]["sample_s"]])
            case_bam = ','.join([str(self.cuffquantpath / x / "abundances.cxb") for x in group_dic["case"]["sample_s"]])

            #cont_bam = ','.join([str(self.tophatpath / x / "accepted_hits.bam") for x in group_dic["cont"]["sample_s"]])
            #case_bam = ','.join([str(self.tophatpath / x / "accepted_hits.bam") for x in group_dic["case"]["sample_s"]])
            _cmd.append(cont_bam)
            _cmd.append(case_bam)
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def generate_md5sum_summay(self, extension):
        outfh = self.md5sum_summary.open("w")
        for file_path in self.readypath.glob(f"**/*.{extension}"):
            if file_path.is_file():
                result = subprocess.run(["md5sum", str(file_path)], capture_output=True, text=True)
                outfh.write(result.stdout)
        outfh.close()
        return 1
    
    def make_bismark_align_sh(self):
        outfh = self.bismark_sh.open("w")
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['bismark']
            _cmd.append("--genome")
            _cmd.append(self.meta_dic["bismark_genome"])
            _cmd.append("-1")
            _cmd.append(str(self.cleanpath / f"{target_id}_R1.fastq.gz"))
            _cmd.append("-2")
            _cmd.append(str(self.cleanpath / f"{target_id}_R2.fastq.gz"))
            _cmd.append("--bowtie2")
            _cmd.append("--nucleotide_coverage")
            _cmd.append("--parallel")
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append("--rg_tag")
            _cmd.append("--rg_id")
            _cmd.append(target_id)
            _cmd.append("--rg_sample")
            _cmd.append(target_id)
            _cmd.append("--output_dir")
            _cmd.append(str(self.bismarkpath))
            _cmd.append("--temp_dir")
            _cmd.append(str(self.bismarkpath / "temp"))
            _cmd.append("--gzip")
            _cmd.append(">")
            _cmd.append(str(self.bismarkpath / f"{target_id}_bismark_align.log"))
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_bismark_dedup_sh(self):
        outfh = self.bismark_sh.open("a")
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['deduplicate_bismark']
            _cmd.append("--paired")
            _cmd.append("--bam")
            _cmd.append("--output_dir")
            _cmd.append(str(self.bismarkpath))
            _cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.bam"))
            _cmd.append(">")
            _cmd.append(str(self.bismarkpath / f"{target_id}_deduplicate_bismark.log"))
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()
    
    def make_bismark_methylcall_sh(self):
        outfh = self.bismark_sh.open("a")
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['bismark_methylation_extractor']
            _cmd.append("--paired-end")
            _cmd.append("--no_overlap")
            _cmd.append("--comprehensive")
            _cmd.append("--CX_context")
            _cmd.append("--gzip")
            _cmd.append("--bedGraph")
            _cmd.append("--cytosine_report")
            _cmd.append("--report")
            _cmd.append("--genome_folder")
            _cmd.append(self.meta_dic["bismark_genome"])
            _cmd.append("--output_dir")
            _cmd.append(str(self.bismarkpath))
            _cmd.append("--multicore")
            _cmd.append(self.meta_dic["cpu"])
            #_cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.bam"))
            _cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.deduplicated.bam"))
            _cmd.append(">")
            _cmd.append(str(self.bismarkpath / f"{target_id}_bismark_methylation_extractor.log"))
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()
            
    def make_bismark_report_sh(self):
        outfh = self.bismark_sh.open("a")
        _cmd = ['pushd']
        _cmd.append(str(self.bismarkpath))
        _cmd.append('&&')
        _cmd.append('bismark2report')
        _cmd.append(">")
        _cmd.append("bismark2report.log")
        _cmd.append('&&')
        _cmd.append('popd')
        outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_bismark_summary_sh(self):
        outfh = self.bismark_sh.open("a")
        _cmd = ['pushd']
        _cmd.append(str(self.bismarkpath))
        _cmd.append('&&')
        _cmd.append('bismark2summary')
        _cmd.append(">")
        _cmd.append("bismark2summary.log")
        _cmd.append('&&')
        _cmd.append('popd')
        outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_bismark_hsmetrics_sh(self):
        outfh = self.bismark_sh.open("a")
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['picard']
            _cmd.append("CollectHsMetrics")
            _cmd.append("--INPUT")
            #_cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.bam"))
            _cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.deduplicated.bam"))
            _cmd.append("--OUTPUT")
            #_cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.bam.hsmetrics.txt"))
            _cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.deduplicated.bam.hsmetrics.txt"))
            _cmd.append("--REFERENCE_SEQUENCE")
            _cmd.append(self.meta_dic["genome_fasta"])
            _cmd.append("--BAIT_INTERVALS")
            _cmd.append(self.meta_dic["interval_list"])
            _cmd.append("--TARGET_INTERVALS")
            _cmd.append(self.meta_dic["interval_list"])
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_bismark_sort_sh(self):
        outfh = self.bismark_sh.open("a")
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['samtools']
            _cmd.append("sort")
            _cmd.append("-@")
            _cmd.append(self.meta_dic["cpu"])
            _cmd.append("-o")
            #_cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.sorted.bam"))
            _cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.deduplicated.sorted.bam"))
            #_cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.bam"))
            _cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.deduplicated.bam"))
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

    def make_bismark_qualimap_sh(self):
        outfh = self.bismark_sh.open("a")
        for target_id, read_dic in self.fastq_dic['ready'].items():
            _cmd = ['qualimap']
            _cmd.append("bamqc")
            _cmd.append("-bam")
            #_cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.sorted.bam"))
            _cmd.append(str(self.bismarkpath / f"{target_id}_R1_bismark_bt2_pe.deduplicated.sorted.bam"))
            _cmd.append("-outdir")
            _cmd.append(str(self.bismarkpath / f"{target_id}_bismark_qualimap"))
            _cmd.append("--paint-chromosome-limits")
            _cmd.append("--feature-file")
            _cmd.append(self.meta_dic["target_bed"])
            _cmd.append("--java-mem-size=4G")
            outfh.write("{0}\n".format(' '.join(_cmd)))
        outfh.close()

def main(args):

    obj = Demux2Ready(args.conf_fn, args.work_dir, args.find_fastq_mode)

    if args.mode in ['mksh_d2r']:
        obj.make_demux2ready_sh()
    elif args.mode in ['prepare']:
        obj.prepare_now()

    elif args.mode in ['mksh_fastp']:
        obj.make_fastp_sh()
    elif args.mode in ['parse_fastp']:
        obj.parse_fastp()

    elif args.mode in ['md5sum']:
        obj.generate_md5sum_summay("fastq.gz")

    elif args.mode in ['mkconf_nfcore_rnaseq']:
        obj.make_nfcore_rnaseq_samplesheet()
    elif args.mode in ['mkconf_nfcore_differentialabundance']:
        obj.make_nfcore_diffabdc_samplesheet()
        obj.make_nfcore_diffabdc_contrasts()
    elif args.mode in ['mkconf_atgcu_rnaseq']:
        obj.make_atgcu_rnaseq_sample()
        obj.make_atgcu_rnaseq_compare()
    elif args.mode in ['samplefile']:
        obj.make_samplefile()

    elif args.mode in ['mksh_cutadapt']: # read length trimming
        obj.cleanpath.mkdir(exist_ok=True)
        obj.make_cutadapt_sh()
    elif args.mode in ['mksh_fastq_screen']:
        obj.fastqscreenpath.mkdir(exist_ok=True)
        obj.make_fastq_screen_sh()
    elif args.mode in ['mksh_fastp_clean']:
        obj.cleanpath.mkdir(exist_ok=True)
        obj.make_fastp_clean_sh()
    elif args.mode in ['mksh_fastp_split']:
        obj.baseadjustpath.mkdir(exist_ok=True)
        obj.splitpath.mkdir(exist_ok=True)
        obj.make_fastp_split_sh()
        obj.make_baseadjust_sh()
    elif args.mode in ['mksh_tophat']:
        obj.tophatpath.mkdir(exist_ok=True)
        obj.make_tophat_sh()
    elif args.mode in ['mksh_cuffquant']:
        obj.cuffquantpath.mkdir(exist_ok=True)
        obj.make_cuffquant_sh()
    elif args.mode in ['mksh_cufflinks']:
        obj.cufflinkspath.mkdir(exist_ok=True)
        obj.make_cufflinks_sh()
    elif args.mode in ['mksh_cuffnorm']:
        obj.cuffnormpath.mkdir(exist_ok=True)
        obj.make_cuffnorm_sh()
    elif args.mode in ['mksh_cuffdiff']:
        obj.cuffdiffpath.mkdir(exist_ok=True)
        obj.make_cuffdiff_sh()
    elif args.mode in ['mksh_bismark']:
        obj.bismarkpath.mkdir(exist_ok=True)
        obj.make_bismark_align_sh() # always first
        obj.make_bismark_dedup_sh()
        obj.make_bismark_methylcall_sh()
        obj.make_bismark_report_sh()
        obj.make_bismark_summary_sh()

        obj.make_bismark_sort_sh()
        obj.make_bismark_hsmetrics_sh()
        obj.make_bismark_qualimap_sh()





if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--conf-fn', default='etc/demux2ready.conf')
    parser.add_argument('--work-dir', default='./')
    parser.add_argument('--find-fastq-mode', choices=('demux','outsourcing'), default='demux')
    parser.add_argument('--mode', choices=('mksh_d2r', 'prepare', 'md5sum',
                                           'mksh_fastp', 'parse_fastp',
                                           'mkconf_nfcore_rnaseq',
                                           'mkconf_nfcore_differentialabundance',
                                           'mkconf_atgcu_rnaseq',
                                           'samplefile',
                                           'mksh_fastp_clean',
                                           'mksh_fastp_split',
                                           'mksh_fastq_screen',
                                           'mksh_cutadapt',
                                           'mksh_tophat', 'mksh_cuffquant',
                                           'mksh_cufflinks', 'mksh_cuffnorm',
                                           'mksh_cuffdiff',
                                           'mksh_bismark'),
                        default='mksh_d2r')
    args = parser.parse_args()
    main(args)

