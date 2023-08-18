

import math
from pathlib import Path
import logging
logging.basicConfig(level=logging.DEBUG)
#logging.basicConfig(level=logging.INFO)

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
        outfh.write("{0}\n".format('\t'.join([str(x) for x in headers])))
        for sample_name, info_dic in self.align_summary_info_dic.items():
            items = [sample_name]
            items.append(info_dic["total"]["input"])
            items.append(info_dic["total"]["mapped"])
            items.append(round((info_dic["total"]["mapped"]/(info_dic["total"]["input"]*1.0))*100, 2))
            items.append(info_dic["total"]["multimap"])
            items.append(round((info_dic["total"]["multimap"]/(info_dic["total"]["input"]*1.0))*100, 2))
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
        self.find_fpkm_tracking_path("genes", "cufflinks")
        self.parse_fpkm_tracking_path("genes")
        if self.gene2desc_fp.is_file():
            self.anno_gene2desc("BIOTYPE", "biotype")
            self.anno_gene2desc("SYMBOL", "symbol")
            self.anno_gene2desc("DESCRIPTION", "description")
        self.write_fpkm_tsv("genes", 0.0, self.wkdir / f"{outprefix}.genes.all.fpkm.tsv")
        self.write_fpkm_tsv("genes", 0.3, self.wkdir / f"{outprefix}.genes.exp.fpkm.tsv")

        # parse isoforms.fpkm_tracking files | cufflinks with -G option
        self.tran_path_dic = dict()
        self.tran_info_dic = dict()
        self.tran_fpkm_dic = dict()
        self.find_fpkm_tracking_path("isoforms", "cufflinks")
        self.parse_fpkm_tracking_path("isoforms")
        self.write_fpkm_tsv("isoforms", 0.0, self.wkdir / f"{outprefix}.trans.all.fpkm.tsv")
        self.write_fpkm_tsv("isoforms", 0.3, self.wkdir / f"{outprefix}.trans.exp.fpkm.tsv")

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

        if molecule_type in ["genes"]:
            self.gene_info_dic = _info_dic
            self.gene_fpkm_dic = _fpkm_dic
        elif molecule_type in ["isoforms"]:
            self.tran_info_dic = _info_dic
            self.tran_fpkm_dic = _fpkm_dic

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

        return 1

class Cuffdiff:
    def __init__(self, wkdir, outprefix, log2fc_criteria, p_criteria):
        self.wkdir = Path(wkdir)
        self.outprefix = outprefix
        self.comp_names = list()
        self._diff_col_s = ['Status', 'CtrlName','CaseName',
                            'CtrlFPKM','CaseFPKM',
                            'Log2FC','p','q', 'UpDown','IsDEG']
        self._info_col_s = ['GeneName', 'Locus']
        self.log2fc_criteria = log2fc_criteria
        self.p_criteria = p_criteria

        # parse gene_exp.diff files
        self.gene_path_dic = dict()
        self.gene_info_dic = dict()
        self.gene_diff_dic = dict()
        self.find_exp_diff_path("gene", "cuffdiff")
        self.parse_exp_diff_path("gene")

        # parse isoform_exp_diff files
        self.tran_path_dic = dict()
        self.tran_info_dic = dict()
        self.tran_diff_dic = dict()
        self.find_exp_diff_path("isoform", "cuffdiff")
        self.parse_exp_diff_path("isoform")

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

    def anno_gene2desc(self, gene2desc_fp, colname, key):
        for line in gene2desc_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ["ID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic["ID"]]
            if _id in self.gene_info_dic:
                self.gene_info_dic[_id].setdefault(key, items[idx_dic[colname]])

        self._info_col_s.append(key)

        return 1

    def write_deg_tsv(self, molecule_type, _outfn):
        if molecule_type in ["gene"]:
            _info_dic = self.gene_info_dic
            _diff_dic = self.gene_diff_dic
        elif molecule_type in ["isoform"]:
            _info_dic = self.tran_info_dic
            _diff_dic = self.tran_diff_dic

        headers = ["GeneID"]
        for comp_name in self.comp_names:
            for _diff_col in self._diff_col_s:
                headers.append('{0}:{1}'.format(comp_name, _diff_col))
        headers.extend(self._info_col_s)

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
            for _info_col in self._info_col_s:
                if _info_col in _info_dic[_id]:
                    items.append(_info_dic[_id][_info_col])
                else:
                    items.append("None")
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1







def main(args):
    logging.info(args)

    if args.biotool in ['tophat']:
        tophat = Tophat(args.wkdir, args.outprefix)

    elif args.biotool in ['cufflinks']:
        cufflinks = Cufflinks(args.wkdir, args.outprefix, args.gene2desc)

    elif args.biotool in ['cuffdiff']:
        cuffdiff = Cuffdiff(args.wkdir, args.outprefix, 0.584, 0.05)
        if Path(args.gene2desc).is_file():
            cuffdiff.anno_gene2desc(Path(args.gene2desc), "SYMBOL", "GeneSymbol")
            cuffdiff.anno_gene2desc(Path(args.gene2desc), "DESCRIPTION", "Description")
            cuffdiff.anno_gene2desc(Path(args.gene2desc), "BIOTYPE", "Biotype")
        cuffdiff.write_deg_tsv("gene", "genes.deg.tsv")
        cuffdiff.write_deg_tsv("isoform", "trans.deg.tsv")



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description="MultiParser")
    parser.add_argument("--wkdir", default="./")
    subparsers = parser.add_subparsers(required=True, dest="biotool")

    subparser = subparsers.add_parser("tophat")
    subparser.add_argument("--outprefix", default="MultiParser.Tophat")

    subparser = subparsers.add_parser("cufflinks")
    subparser.add_argument("--outprefix", default="MultiParser.Cufflinks")
    subparser.add_argument("--gene2desc", default="anno.gene2desc.tsv")

    subparser = subparsers.add_parser("cuffdiff")
    subparser.add_argument("--outprefix", default="MultiParser.Cuffdiff")
    subparser.add_argument("--gene2desc", default="anno.gene2desc.tsv")

    args = parser.parse_args()
    main(args)






