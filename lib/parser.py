

import math
import json
from pathlib import Path
import logging
logging.basicConfig(level=logging.DEBUG)
#logging.basicConfig(level=logging.INFO)


class Fastp:
    def __init__(self, wkdir, outprefix):
        self.wkdir = Path(wkdir)
        self.outprefix = outprefix
        self.json_path_dic = dict()

    def find_json_path(self, target_dir):
        for json_path in Path(target_dir).glob("*.fastp.json"):
            sample_name = json_path.parts[-1].replace(".fastp.json", "")
            logging.info(f"Found {str(json_path)} of {sample_name}")
            self.json_path_dic.setdefault(sample_name, json.load(json_path.open("r")))
        return 1

    def write_summary_tsv(self, tag):
        outfp = self.wkdir / f"{self.outprefix}.tsv"
        logging.info(f"writing fastp summary table... {str(outfp)}")
        outfh = outfp.open("w")
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
        for sample_name, info_dic in self.json_path_dic.items():
            items = [sample_name]
            total_reads = info_dic["summary"][tag]["total_reads"]
            items.append(total_reads)
            items.append(info_dic["summary"][tag]["total_bases"])
            items.append(round(int(info_dic["summary"][tag]["total_bases"])*0.000000001, 3))
            items.append(round(float(info_dic["summary"][tag]["gc_content"])*100, 3))
            items.append(round(float(info_dic["summary"][tag]["q20_rate"])*100, 3))
            items.append(round(float(info_dic["summary"][tag]["q30_rate"])*100, 3))
            items.append(round(float(info_dic["duplication"]["rate"])*100, 3))
            if "insert_size" in info_dic:
                items.append(info_dic["insert_size"]["peak"])
            else:
                items.append("-")
            items.append(info_dic["summary"][tag]["read1_mean_length"])
            if "read2_mean_length" in info_dic["summary"][tag]:
                items.append(info_dic["summary"][tag]["read2_mean_length"])
            else:
                items.append("-")
            items.append(round(int(info_dic["filtering_result"]["passed_filter_reads"])/float(total_reads)*100, 3))
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1


class Tophat:
    def __init__(self, wkdir, outprefix):
        self.wkdir = Path(wkdir)

        # parse align_summary.txt files
        self.align_summary_path_dic = dict()
        self.align_summary_info_dic = dict()
        self.find_align_summary_path("tophat")
        self.parse_align_summary_path()
        self.write_align_summary_tsv(self.wkdir / f"{outprefix}.align_summary.tsv")

    def write_align_summary_tsv(self, outfp):
        logging.info(f"writing align summary table... {str(outfp)}")
        outfh = outfp.open("w")
        headers = ["SampleName"]
        headers.append("InputCount")
        headers.append("MappedCount")
        headers.append("MappedRate")
        headers.append("MultiMapCount")
        headers.append("MultiMapRate")
        headers.append("AlignedAsPairCount")
        headers.append("AlignedAsPairRate")
        headers.append("AlignedAsPair-MultiRate")
        headers.append("AlignedAsPair-DiscordantRate")
        headers.append("AlignedAsPair-ConcordantRate")
        outfh.write("{0}\n".format('\t'.join([str(x) for x in headers])))
        for sample_name, info_dic in self.align_summary_info_dic.items():
            items = [sample_name]
            items.append(info_dic["total"]["input"])
            items.append(info_dic["total"]["mapped"])
            items.append(round((info_dic["total"]["mapped"]/(info_dic["total"]["input"]*1.0))*100, 2))
            items.append(info_dic["total"]["multimap"])
            items.append(round((info_dic["total"]["multimap"]/(info_dic["total"]["mapped"]*1.0))*100, 2))
            items.append(info_dic["total"]["aligned-pair"])
            items.append(round((info_dic["total"]["aligned-pair"]/(info_dic["total"]["mapped"]*1.0))*100, 2))
            items.append(round((info_dic["total"]["aligned-pair-multi"]/(info_dic["total"]["aligned-pair"]*1.0))*100, 2))
            items.append(round((info_dic["total"]["aligned-pair-discordant"]/(info_dic["total"]["aligned-pair"]*1.0))*100, 2))
            items.append(round(info_dic["total"]["aligned-pair-concordant-rate"], 2))
            logging.info(items)
            outfh.write("{0}\n".format('\t'.join([str(x) for x in items])))
        outfh.close()
        return 1


    def parse_align_summary_path(self):
        for sample_name, align_summary_path in self.align_summary_path_dic.items():
            logging.debug(sample_name)
            metrix = list()
            for line in align_summary_path.open():
                items = line.strip().split()
                if not items:
                    continue
                metrix.append(items)
            logging.debug(len(metrix))
            for items in metrix:
                logging.debug(items)

            _dic = dict()
            _dic.setdefault("r1", {}).setdefault("input", 0)
            _dic.setdefault("r1", {}).setdefault("mapped", 0)
            _dic.setdefault("r1", {}).setdefault("multimap", 0)
            _dic.setdefault("r2", {}).setdefault("input", 0)
            _dic.setdefault("r2", {}).setdefault("mapped", 0)
            _dic.setdefault("r2", {}).setdefault("multimap", 0)
            _dic.setdefault("total", {}).setdefault("input", 0)
            _dic.setdefault("total", {}).setdefault("mapped", 0)
            _dic.setdefault("total", {}).setdefault("multimap", 0)
            _dic.setdefault("total", {}).setdefault("aligned-pair", 0)
            _dic.setdefault("total", {}).setdefault("aligned-pair-multi", 0)
            _dic.setdefault("total", {}).setdefault("aligned-pair-discordant", 0)
            _dic.setdefault("total", {}).setdefault("aligned-pair-concordant-rate", 0.0)
            if len(metrix) in [13]:
                _dic["r1"]["input"] = int(metrix[1][2])
                _dic["r1"]["mapped"] = int(metrix[2][2])
                _dic["r1"]["multimap"] = int(metrix[3][2])
                _dic["r2"]["input"] = int(metrix[5][2])
                _dic["r2"]["mapped"] = int(metrix[6][2])
                _dic["r2"]["multimap"] = int(metrix[7][2])
                _dic["total"]["input"] = _dic["r1"]["input"] + _dic["r2"]["input"]
                _dic["total"]["mapped"] = _dic["r1"]["mapped"] + _dic["r2"]["mapped"]
                _dic["total"]["multimap"] = _dic["r1"]["multimap"] + _dic["r2"]["multimap"]
                _dic["total"]["aligned-pair"] = int(metrix[9][2])*2
                _dic["total"]["aligned-pair-multi"] = int(metrix[10][2])*2
                _dic["total"]["aligned-pair-discordant"] = int(metrix[11][0])*2
                _dic["total"]["aligned-pair-concordant-rate"] = float(metrix[12][0].replace("%", ""))
            else:
                logging.critical("It is need to dev for this.")

            self.align_summary_info_dic.setdefault(sample_name, _dic)
        logging.debug(self.align_summary_info_dic)
        return 1

    def find_align_summary_path(self, tophat_dir):
        for align_summary_path in Path(tophat_dir).glob("*/align_summary.txt"):
            sample_name = align_summary_path.parts[-2]
            logging.info(f"Found {str(align_summary_path)} of {sample_name}")
            self.align_summary_path_dic.setdefault(sample_name, align_summary_path)
        return 1


class Cufflinks:
    def __init__(self, wkdir, outprefix, gene2desc):
        self.wkdir = Path(wkdir)
        self.sample_names = list()
        self.gene2desc_fp = Path(gene2desc)

        # setting vars for annotation
        self.header_names = list()
        self.anno_names = list()
        _names = [("locus", "Locus"),
                  ("name", "GeneName"),
                  ("symbol", "GeneSymbol"),
                  ("description", "Description"),
                  ("biotype", "Biotype")]
        for anno_name, header_name in _names:
            self.anno_names.append(anno_name)
            self.header_names.append(header_name)

        # parse gene.fpkm_tracking files | cufflinks with -G option
        self.gene_path_dic = dict()
        self.gene_info_dic = dict()
        self.gene_fpkm_dic = dict()
        self.gene_stat_dic = dict()
        self.find_fpkm_tracking_path("genes", "cufflinks")
        self.parse_fpkm_tracking_path("genes")
        if self.gene2desc_fp.is_file():
            self.anno_gene2desc("BIOTYPE", "biotype")
            self.anno_gene2desc("SYMBOL", "symbol")
            self.anno_gene2desc("DESCRIPTION", "description")
        self.write_fpkm_tsv("genes", 0.0, self.wkdir / f"{outprefix}.genes.all.fpkm.tsv")
        self.write_fpkm_tsv("genes", 0.3, self.wkdir / f"{outprefix}.genes.exp.fpkm.tsv")
        self.write_stat_exp("genes", self.wkdir / f"{outprefix}.genes.stat.expressed.tsv")
        self.write_stat_coexp("genes", self.wkdir / f"{outprefix}.genes.stat.co-expressed.tsv")

        # parse isoforms.fpkm_tracking files | cufflinks with -G option
        self.tran_path_dic = dict()
        self.tran_info_dic = dict()
        self.tran_fpkm_dic = dict()
        self.tran_stat_dic = dict()
        self.find_fpkm_tracking_path("isoforms", "cufflinks")
        self.parse_fpkm_tracking_path("isoforms")
        self.write_fpkm_tsv("isoforms", 0.0, self.wkdir / f"{outprefix}.trans.all.fpkm.tsv")
        self.write_fpkm_tsv("isoforms", 0.3, self.wkdir / f"{outprefix}.trans.exp.fpkm.tsv")
        self.write_stat_exp("isoforms", self.wkdir / f"{outprefix}.trans.stat.expressed.tsv")
        self.write_stat_coexp("isoforms", self.wkdir / f"{outprefix}.trans.stat.co-expressed.tsv")

    def write_stat_coexp(self, molecule_type, outfp):
        if molecule_type in ["genes"]:
            _stat_dic = self.gene_stat_dic
            _fpkm_dic = self.gene_fpkm_dic
        elif molecule_type in ["isoforms"]:
            _stat_dic = self.tran_stat_dic
            _fpkm_dic = self.tran_fpkm_dic

        logging.info(f"writing STAT of Co-Expressed gene table... {str(outfp)}")

        _stat_col_s = ["fpkm>=0.3", "fpkm>=1", "fpkm>=10", "fpkm>=100", "fpkm>=1000"]

        _coexp_dic = dict()
        for i in range(len(self.sample_names)+1):
            for _stat_col in _stat_col_s:
                _coexp_dic.setdefault(i, {}).setdefault(_stat_col, 0)

        for _id, info_dic in _fpkm_dic.items():
            _fpkm_s = list()
            for sample_name in self.sample_names:
                _fpkm_s.append(float(info_dic[sample_name]))
            for _stat_col in _stat_col_s:
                _fpkm_cutoff = 0.0
                if _stat_col in ["fpkm>=0.3"]:
                    _fpkm_cutoff = 0.3
                elif _stat_col in ["fpkm>=1"]:
                    _fpkm_cutoff = 1
                elif _stat_col in ["fpkm>=10"]:
                    _fpkm_cutoff = 10
                elif _stat_col in ["fpkm>=100"]:
                    _fpkm_cutoff = 100
                elif _stat_col in ["fpkm>=1000"]:
                    _fpkm_cutoff = 1000
                _isexp_s = ["O" if x >= _fpkm_cutoff else "X" for x in _fpkm_s]
                _coexp_dic[_isexp_s.count("O")][_stat_col] += 1

        headers = ["CoExpSampleCount"]
        headers.extend(_stat_col_s)
        outfh = outfp.open("w")
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for coexp_sample_count, info_dic in sorted(_coexp_dic.items(), reverse=True):
            items = [coexp_sample_count]
            for _stat_col in _stat_col_s:
                items.append(info_dic[_stat_col])
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1

    def write_stat_exp(self, molecule_type, outfp):
        if molecule_type in ["genes"]:
            _stat_dic = self.gene_stat_dic
        elif molecule_type in ["isoforms"]:
            _stat_dic = self.tran_stat_dic

        logging.info(f"writing STAT of Expressed gene table... {str(outfp)}")

        _stat_col_s = ["fpkm>=0.3", "fpkm>=1", "fpkm>=10",
                       "fpkm>=100", "fpkm>=1000"]
        headers = ["SampleNames"]
        headers.extend(_stat_col_s)
        outfh = outfp.open("w")
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for sample_name in self.sample_names:
            info_dic = _stat_dic[sample_name]
            items = [sample_name]
            for _stat_col in _stat_col_s:
                items.append(len(info_dic[_stat_col]))
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1


    def write_fpkm_tsv(self, molecule_type, fpkm_criteria, outfp):
        if molecule_type in ["genes"]:
            _info_dic = self.gene_info_dic
            _fpkm_dic = self.gene_fpkm_dic
        elif molecule_type in ["isoforms"]:
            _info_dic = self.tran_info_dic
            _fpkm_dic = self.tran_fpkm_dic

        logging.info(f"writing FPKM table... {str(outfp)}")

        headers = ["GeneID"]
        headers.extend(self.sample_names)
        headers.extend(self.header_names)

        outfh = outfp.open("w")
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for _id, _sub_dic in _fpkm_dic.items():
            items = [_id]
            fpkm_s = list()
            for sample_name in self.sample_names:
                _fpkm = _sub_dic[sample_name]
                fpkm_s.append(_fpkm)
                items.append(round(_fpkm, 4))
                #items.append(int(_fpkm))

            for anno_name in self.anno_names:
                if anno_name in _info_dic[_id]:
                    items.append(_info_dic[_id][anno_name])
                else:
                    items.append("None")

            if max(fpkm_s) >= fpkm_criteria:
                outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1

    def anno_gene2desc(self, colname, key):
        for line in self.gene2desc_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ["ID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic["ID"]]
            if _id in self.gene_info_dic:
                self.gene_info_dic[_id].setdefault(key, items[idx_dic[colname]])
        return 1

    def parse_fpkm_tracking_path(self, molecule_type):
        if molecule_type in ["genes"]:
            _path_dic = self.gene_path_dic
        elif molecule_type in ["isoforms"]:
            _path_dic = self.tran_path_dic

        _info_dic = dict()
        _fpkm_dic = dict()
        _stat_dic = dict()
        for sample_name, fpkm_tracking_path in _path_dic.items():
            for line in fpkm_tracking_path.open():
                items = line.rstrip('\n').split('\t')
                if items[0] in ["tracking_id"]:
                    idx_dic = dict()
                    for idx, item in enumerate(items):
                        idx_dic.setdefault(item, idx)
                    continue
                _id = items[idx_dic["tracking_id"]]
                _name = items[idx_dic["gene_short_name"]]
                _locus = items[idx_dic["locus"]]
                _fpkm = float(items[idx_dic["FPKM"]])

                _info_dic.setdefault(_id, {}).setdefault("name", _name)
                _info_dic.setdefault(_id, {}).setdefault("locus", _locus)
                _fpkm_dic.setdefault(_id, {}).setdefault(sample_name, _fpkm)
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=0.3", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=1", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=10", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=100", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=1000", [])

                if float(_fpkm) >= 1000:
                    _stat_dic[sample_name]["fpkm>=1000"].append(_id)
                if float(_fpkm) >= 100:
                    _stat_dic[sample_name]["fpkm>=100"].append(_id)
                if float(_fpkm) >= 10:
                    _stat_dic[sample_name]["fpkm>=10"].append(_id)
                if float(_fpkm) >= 1.0:
                    _stat_dic[sample_name]["fpkm>=1"].append(_id)
                if float(_fpkm) >= 0.3:
                    _stat_dic[sample_name]["fpkm>=0.3"].append(_id)


        if molecule_type in ["genes"]:
            self.gene_info_dic = _info_dic
            self.gene_fpkm_dic = _fpkm_dic
            self.gene_stat_dic = _stat_dic
        elif molecule_type in ["isoforms"]:
            self.tran_info_dic = _info_dic
            self.tran_fpkm_dic = _fpkm_dic
            self.tran_stat_dic = _stat_dic

        return 1

    def find_fpkm_tracking_path(self, molecule_type, cufflinks_dir):
        for fpkm_tracking_path in Path(cufflinks_dir).glob(f"*/{molecule_type}.fpkm_tracking"):
            sample_name = fpkm_tracking_path.parts[-2]
            if sample_name not in self.sample_names:
                self.sample_names.append(sample_name)
            logging.info(f"Found {str(fpkm_tracking_path)} of {sample_name}")

            if molecule_type in ["genes"]:
                self.gene_path_dic.setdefault(sample_name, fpkm_tracking_path)
            elif molecule_type in ["isoforms"]:
                self.tran_path_dic.setdefault(sample_name, fpkm_tracking_path)
        self.sample_names.sort()

        return 1

class Cuffdiff:
    def __init__(self, wkdir, outprefix, log2fc_criteria, p_criteria):
        self.wkdir = Path(wkdir)
        self.outprefix = outprefix
        self.comp_names = list()
        self._diff_col_s = ['Status', 'CtrlName','CaseName',
                            'CtrlFPKM','CaseFPKM',
                            'Log2FC','p','q', 'UpDown','IsDEG']
        self.gene_info_col_s = ['GeneName', 'Locus']
        self.tran_info_col_s = ['GeneName', 'Locus']
        self.log2fc_criteria = log2fc_criteria
        self.p_criteria = p_criteria

        # parse gene_exp.diff files
        self.gene_path_dic = dict()
        self.gene_info_dic = dict()
        self.gene_diff_dic = dict()
        self.gene_stat_dic = dict()
        self.find_exp_diff_path("gene", "cuffdiff")
        self.parse_exp_diff_path("gene")
        self.stat_exp_diff("gene")

        # parse isoform_exp_diff files
        self.tran_path_dic = dict()
        self.tran_info_dic = dict()
        self.tran_diff_dic = dict()
        self.find_exp_diff_path("isoform", "cuffdiff")
        self.parse_exp_diff_path("isoform")
        self.stat_exp_diff("isoform")

    def find_exp_diff_path(self, molecule_type, cuffdiff_dir):
        for exp_diff_path in Path(cuffdiff_dir).glob(f"*/{molecule_type}_exp.diff"):
            comp_name = exp_diff_path.parts[-2]
            if comp_name not in self.comp_names:
                self.comp_names.append(comp_name)
            logging.info(f"Found {str(exp_diff_path)} of {comp_name}")

            if molecule_type in ["gene"]:
                self.gene_path_dic.setdefault(comp_name, exp_diff_path)
            elif molecule_type in ["isoform"]:
                self.tran_path_dic.setdefault(comp_name, exp_diff_path)

        self.comp_names.sort()

        return 1

    def stat_exp_diff(self, molecule_type):
        if molecule_type in ["gene"]:
            _diff_dic = self.gene_diff_dic
        elif molecule_type in ["isoform"]:
            _diff_dic = self.tran_diff_dic

        _stat_dic = dict()
        for _id, _sub_diff_dic in _diff_dic.items():
            for comp_name, info_dic in _sub_diff_dic.items():

                _stat_dic.setdefault(comp_name, {}).setdefault('DEG', 0)
                _stat_dic.setdefault(comp_name, {}).setdefault('DEG_Up', 0)
                _stat_dic.setdefault(comp_name, {}).setdefault('DEG_Down', 0)
                if info_dic['UpDown'] in ['up'] and info_dic['IsDEG'] in ['yes']:
                    _stat_dic[comp_name]['DEG'] += 1
                    _stat_dic[comp_name]['DEG_Up'] += 1
                elif info_dic['UpDown'] in ['down'] and info_dic['IsDEG'] in ['yes']:
                    _stat_dic[comp_name]['DEG'] += 1
                    _stat_dic[comp_name]['DEG_Down'] += 1

                _stat_dic.setdefault(comp_name, {}).setdefault('Delta_Pstv', 0)
                _stat_dic.setdefault(comp_name, {}).setdefault('Delta_Ngtv', 0)
                if info_dic['UpDown'] in ['up']:
                    _stat_dic[comp_name]['Delta_Pstv'] += 1
                elif info_dic['UpDown'] in ['down']:
                    _stat_dic[comp_name]['Delta_Ngtv'] += 1

        if molecule_type in ["gene"]:
            self.gene_stat_dic = _stat_dic
        elif molecule_type in ["isoform"]:
            self.tran_stat_dic = _stat_dic
        return 1


    def parse_exp_diff_path(self, molecule_type):
        if molecule_type in ["gene"]:
            _path_dic = self.gene_path_dic
        elif molecule_type in ["isoform"]:
            _path_dic = self.tran_path_dic

        _info_dic = dict()
        _diff_dic = dict()
        for comp_name, exp_diff_path in _path_dic.items():
            for line in exp_diff_path.open():
                items = line.rstrip('\n').split('\t')
                if items[0] in ["test_id"]:
                    idx_dic = dict()
                    for idx, item in enumerate(items):
                        idx_dic.setdefault(item, idx)
                    continue
                _id = items[idx_dic["test_id"]]
                _name = items[idx_dic["gene"]]
                _locus = items[idx_dic["locus"]]
                ctrl_id = items[idx_dic['sample_1']]
                case_id = items[idx_dic['sample_2']]
                status = items[idx_dic['status']]
                ctrl_value = items[idx_dic['value_1']]
                case_value = items[idx_dic['value_2']]
                log2fc = math.log(float(case_value)+1.0,2) - math.log(float(ctrl_value)+1.0,2)
                p = items[idx_dic['p_value']]
                q = items[idx_dic['q_value']]

                _info_dic.setdefault(_id, {}).setdefault("GeneName", _name)
                _info_dic.setdefault(_id, {}).setdefault("Locus", _locus)

                _diff_dic.setdefault(_id, {}).setdefault(comp_name, {})
                _diff_dic[_id][comp_name].setdefault('Status', status)
                _diff_dic[_id][comp_name].setdefault('CtrlName', ctrl_id)
                _diff_dic[_id][comp_name].setdefault('CaseName', case_id)
                _diff_dic[_id][comp_name].setdefault('CtrlFPKM', ctrl_value)
                _diff_dic[_id][comp_name].setdefault('CaseFPKM', case_value)
                _diff_dic[_id][comp_name].setdefault('Log2FC', log2fc)
                if log2fc > 0.0:
                    _diff_dic[_id][comp_name].setdefault('UpDown', 'up')
                elif log2fc < 0.0:
                    _diff_dic[_id][comp_name].setdefault('UpDown', 'down')
                else:
                    _diff_dic[_id][comp_name].setdefault('UpDown', 'flat')
                _diff_dic[_id][comp_name].setdefault('p', p)
                _diff_dic[_id][comp_name].setdefault('q', q)
                if abs(log2fc) >= self.log2fc_criteria and float(p) < self.p_criteria:
                    _diff_dic[_id][comp_name].setdefault('IsDEG', 'yes')
                else:
                    _diff_dic[_id][comp_name].setdefault('IsDEG', 'no')

        if molecule_type in ["gene"]:
            self.gene_info_dic = _info_dic
            self.gene_diff_dic = _diff_dic
        elif molecule_type in ["isoform"]:
            self.tran_info_dic = _info_dic
            self.tran_diff_dic = _diff_dic

        return 1

    def anno_xxxx2desc(self, molecule_type, xxxx2desc_fp, colname, key):
        for line in xxxx2desc_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ["ID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic["ID"]]
            if molecule_type in ["gene"] and _id in self.gene_info_dic:
                self.gene_info_dic[_id].setdefault(key, items[idx_dic[colname]])
            if molecule_type in ["isoform"] and _id in self.tran_info_dic:
                self.tran_info_dic[_id].setdefault(key, items[idx_dic[colname]])

        if molecule_type in ["gene"] and key not in self.gene_info_col_s:
            self.gene_info_col_s.append(key)
        if molecule_type in ["isoform"] and key not in self.tran_info_col_s:
            self.tran_info_col_s.append(key)

        return 1

    def write_deg_tsv(self, molecule_type, _outfn):
        if molecule_type in ["gene"]:
            _info_dic = self.gene_info_dic
            _diff_dic = self.gene_diff_dic
            _info_col_s = self.gene_info_col_s
        elif molecule_type in ["isoform"]:
            _info_dic = self.tran_info_dic
            _diff_dic = self.tran_diff_dic
            _info_col_s = self.tran_info_col_s

        headers = ["GeneID"]
        for comp_name in self.comp_names:
            for _diff_col in self._diff_col_s:
                headers.append('{0}:{1}'.format(comp_name, _diff_col))
        headers.extend(_info_col_s)

        outfp = self.wkdir / f"{self.outprefix}.{_outfn}"
        logging.info(f"writing DEG table... {str(outfp)}")
        outfh = outfp.open("w")
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for _id, _sub_dic in _diff_dic.items():
            items = [_id]
            for comp_name in self.comp_names:
                for _diff_col in self._diff_col_s:
                    if _diff_col in ["CtrlFPKM", "CaseFPKM", "Log2FC"]:
                        items.append(round(float(_sub_dic[comp_name][_diff_col]), 4))
                    else:
                        items.append(_sub_dic[comp_name][_diff_col])
            for _info_col in _info_col_s:
                if _info_col in _info_dic[_id]:
                    items.append(_info_dic[_id][_info_col])
                else:
                    items.append("None")
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1

    def write_deg_stat_tsv(self, molecule_type, _outfn):
        if molecule_type in ["gene"]:
            _stat_dic = self.gene_stat_dic
        elif molecule_type in ["isoform"]:
            _stat_dic = self.tran_stat_dic

        headers = ["ComparisonID"]
        headers.append("DEG")
        headers.append("DEG_Up")
        headers.append("DEG_Down")
        headers.append("Delta_Pstv")
        headers.append("Delta_Ngtv")

        outfp = self.wkdir / f"{self.outprefix}.{_outfn}"
        logging.info(f"writing DEG Stat table... {str(outfp)}")
        outfh = outfp.open("w")
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for comp_name in self.comp_names:
            info_dic = _stat_dic[comp_name]
            items = [comp_name]
            items.append(info_dic["DEG"])
            items.append(info_dic["DEG_Up"])
            items.append(info_dic["DEG_Down"])
            items.append(info_dic["Delta_Pstv"])
            items.append(info_dic["Delta_Ngtv"])
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()

        return 1






def main(args):
    logging.info(args)

    if args.biotool in ['fastp']:

        fastp = Fastp(args.wkdir, args.outprefix)
        fastp.find_json_path(args.target_dir)
        if args.target_dir in ['cleanfastq']:
            fastp.write_summary_tsv("after_filtering")
            #fastp.write_summary_tsv("before_filtering")
        else:
            fastp.write_summary_tsv("before_filtering")


    elif args.biotool in ['tophat']:

        tophat = Tophat(args.wkdir, args.outprefix)

    elif args.biotool in ['cufflinks']:

        cufflinks = Cufflinks(args.wkdir, args.outprefix, args.gene2desc)

    elif args.biotool in ['cuffdiff']:

        cuffdiff = Cuffdiff(args.wkdir, args.outprefix, args.log2fc_criteria, args.p_criteria)

        if Path(args.gene2desc).is_file():
            cuffdiff.anno_xxxx2desc("gene", Path(args.gene2desc), "SYMBOL", "GeneSymbol")
            cuffdiff.anno_xxxx2desc("gene", Path(args.gene2desc), "DESCRIPTION", "Description")
            cuffdiff.anno_xxxx2desc("gene", Path(args.gene2desc), "BIOTYPE", "Biotype")

        if Path(args.tran2desc).is_file():
            cuffdiff.anno_xxxx2desc("isoform", Path(args.tran2desc), "SYMBOL", "TranscriptSymbol")
            cuffdiff.anno_xxxx2desc("isoform", Path(args.tran2desc), "DESCRIPTION", "Description")
            cuffdiff.anno_xxxx2desc("isoform", Path(args.tran2desc), "BIOTYPE", "Biotype")

        cuffdiff.write_deg_tsv("gene", "genes.deg.tsv")
        cuffdiff.write_deg_tsv("isoform", "trans.deg.tsv")
        cuffdiff.write_deg_stat_tsv("gene", "genes.deg.stat.tsv")
        cuffdiff.write_deg_stat_tsv("isoform", "trans.deg.stat.tsv")



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description="MultiParser")
    parser.add_argument("--wkdir", default="./")
    subparsers = parser.add_subparsers(required=True, dest="biotool")

    subparser = subparsers.add_parser("fastp")
    subparser.add_argument("--outprefix", default="MultiParser.Fastp.Raw") # MultiParser.Fastp.Clean if cleanfastq
    subparser.add_argument("--target-dir", default="readyfastq") # cleanfastq

    subparser = subparsers.add_parser("tophat")
    subparser.add_argument("--outprefix", default="MultiParser.Tophat")

    subparser = subparsers.add_parser("cufflinks")
    subparser.add_argument("--outprefix", default="MultiParser.Cufflinks")
    subparser.add_argument("--gene2desc", default="anno.gene2desc.tsv")

    subparser = subparsers.add_parser("cuffdiff")
    subparser.add_argument("--outprefix", default="MultiParser.Cuffdiff")
    subparser.add_argument("--gene2desc", default="anno.gene2desc.tsv")
    subparser.add_argument("--tran2desc", default="anno.transcript2desc.tsv")
    subparser.add_argument("--log2fc-criteria", default=0.584, type=float)
    subparser.add_argument("--p-criteria", default=0.05, type=float)

    args = parser.parse_args()
    main(args)






