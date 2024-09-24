

from pathlib import Path
import gzip
import logging
logging.basicConfig(level=logging.INFO)

class Bismark:
    def __init__(self, wkdir, sample_name, cytosine_context):
        self.wkdir = Path(wkdir)
        self.sample_name = sample_name

        self.cx_report_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.CX_report.txt.gz"
        self.cpg_report_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz"
        self.chg_report_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.CHG_report.txt.gz"
        self.chh_report_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.CHH_report.txt.gz"

        self.bedgraph_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.bedGraph.gz"
        self.cpg_bedgraph_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.CpG.bedGraph"
        self.chg_bedgraph_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.CHG.bedGraph"
        self.chh_bedgraph_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.CHH.bedGraph"

        self.cov_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
        self.cpg_cov_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.CpG.bismark.cov.gz"
        self.chg_cov_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.CHG.bismark.cov.gz"
        self.chh_cov_fp = self.wkdir / f"{self.sample_name}_R1_bismark_bt2_pe.deduplicated.CHH.bismark.cov.gz"

        self.cytosine_context = self.load_cytosine_context(cytosine_context)
        self.pos_dict = {}



    def load_cytosine_context(self, cytosine_context):
        if cytosine_context == "CpG":
            return "CpG"
        elif cytosine_context == "CHG":
            return "CHG"
        elif cytosine_context == "CHH":
            return "CHH"
        else:
            raise ValueError(f"Invalid cytosine context: {cytosine_context} (must be CpG, CHG, or CHH)")

    def split_cx_report(self):
        logging.info("Splitting CpG, CHG, CHH reports")
        outfh_cpg = gzip.open(self.cpg_report_fp, 'wb')
        outfh_chg = gzip.open(self.chg_report_fp, 'wb')
        outfh_chh = gzip.open(self.chh_report_fp, 'wb')
        for line in gzip.open(self.cx_report_fp, 'rt'):
            items = line.rstrip().split("\t")
            if not len(items) in [7]:
                continue
            if items[5] == "CG":
                outfh_cpg.write(line.encode('utf-8'))
            elif items[5] == "CHG":
                outfh_chg.write(line.encode('utf-8'))
            elif items[5] == "CHH":
                outfh_chh.write(line.encode('utf-8'))
        outfh_cpg.close()
        outfh_chg.close()
        outfh_chh.close()
    
    def load_pos_from_cx_report(self):
        logging.info(f"Loading positions from {self.cytosine_context}_report")
        if self.cytosine_context == "CpG":
            cx_report_fp = self.cpg_report_fp
        elif self.cytosine_context == "CHG":
            cx_report_fp = self.chg_report_fp
        elif self.cytosine_context == "CHH":
            cx_report_fp = self.chh_report_fp

        for line in gzip.open(cx_report_fp, 'rt'):
            items = line.rstrip().split("\t")
            if not len(items) in [7]:
                continue
            if int(items[3]) or int(items[4]): # load covered position only
                self.pos_dict.setdefault(str(items[0]), []).append(int(items[1]))

        for chrom, pos_list in self.pos_dict.items():
            logging.info(f"Loaded {chrom}: {len(pos_list)} positions")

    def extract_cx_from_bedgraph(self):
        logging.info(f"Extracted {self.cytosine_context} from bedGraph")
        if not self.pos_dict:
            raise ValueError("Positions not loaded. Call load_pos_from_cx_report first.")

        if self.cytosine_context == "CpG":
            bedgraph_fp = self.cpg_bedgraph_fp
        elif self.cytosine_context == "CHG":
            bedgraph_fp = self.chg_bedgraph_fp
        elif self.cytosine_context == "CHH":
            bedgraph_fp = self.chh_bedgraph_fp

        outfh = open(bedgraph_fp, 'w')
        for line in gzip.open(self.bedgraph_fp, 'rt'):
            if line.startswith("track"):
                #outfh.write(line.encode('utf-8'))
                outfh.write(line)
                continue
            items = line.rstrip().split("\t")
            if int(items[2]) in self.pos_dict[items[0]]:
                #outfh.write(line.encode('utf-8'))
                outfh.write(line)
        outfh.close()

    def extract_cx_from_bismark_cov(self):
        logging.info(f"Extracted {self.cytosine_context} from bismark.cov")
        if not self.pos_dict:
            raise ValueError("Positions not loaded. Call load_pos_from_cx_report first.")

        if self.cytosine_context == "CpG":
            cov_fp = self.cpg_cov_fp
        elif self.cytosine_context == "CHG":
            cov_fp = self.chg_cov_fp
        elif self.cytosine_context == "CHH":
            cov_fp = self.chh_cov_fp

        outfh = gzip.open(cov_fp, 'wb')
        for line in gzip.open(self.cov_fp, 'rt'):
            items = line.rstrip().split("\t")
            if int(items[3]) in self.pos_dict[items[0]]:
                outfh.write(line.encode('utf-8'))
        outfh.close()




def main(args):

    bismark = Bismark(args.wkdir, args.sample_name, args.cytosine_context)
    #bismark.split_cx_report()
    bismark.load_pos_from_cx_report()
    bismark.extract_cx_from_bedgraph()
    #bismark.extract_cx_from_bismark_cov()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--wkdir", type=str, required=True)
    parser.add_argument("--sample_name", type=str, required=True)
    parser.add_argument("--cytosine_context", type=str, required=True,
                        choices=["CpG", "CHG", "CHH"])
    args = parser.parse_args()
    main(args)
