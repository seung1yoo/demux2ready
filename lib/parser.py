

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
    def __init__(self, wkdir, outprefix):
        self.wkdir = Path(wkdir)
        self.sample_names = list()

        # parse gene.fpkm_tracking files | cufflinks with -G option
        self.gene_path_dic = dict()
        self.gene_info_dic = dict()
        self.gene_fpkm_dic = dict()
        self.find_fpkm_tracking_path("genes", "cufflinks")
        self.parse_fpkm_tracking_path("genes")
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
        outfh = outfp.open("w")
        headers = ["GeneID"]
        headers.append("Locus")
        headers.append("GeneName")
        headers.extend(self.sample_names)
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for _id, _sub_dic in _fpkm_dic.items():
            items = [_id]
            items.append(_info_dic[_id]["locus"])
            items.append(_info_dic[_id]["name"])
            fpkm_s = list()
            for sample_name in self.sample_names:
                _fpkm = _sub_dic[sample_name]
                fpkm_s.append(_fpkm)
                items.append(round(_fpkm, 4))
            if max(fpkm_s) >= fpkm_criteria:
                outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
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

def main(args):
    logging.info(args)

    if args.biotool in ['tophat']:
        tophat = Tophat(args.wkdir, args.outprefix)
    elif args.biotool in ['cufflinks']:
        tophat = Cufflinks(args.wkdir, args.outprefix)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description="MultiParser")
    parser.add_argument("--wkdir", default="./")
    subparsers = parser.add_subparsers(required=True, dest="biotool")

    subparser = subparsers.add_parser("tophat")
    subparser.add_argument("--outprefix", default="MultiParser.Tophat")

    subparser = subparsers.add_parser("cufflinks")
    subparser.add_argument("--outprefix", default="MultiParser.Cufflinks")

    args = parser.parse_args()
    main(args)

