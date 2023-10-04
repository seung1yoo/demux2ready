


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import os


class PairedEnd2SingleEnd:

    def __init__(self, outprefix):
        self.r2rc_fn = f"{outprefix}.r2rc.temp.gz"
        self.concat_fn = f"{outprefix}_concat_R1.fastq.gz"

    def do_r2rc(self, in2):
        out_fh = gzip.open(self.r2rc_fn, "wb")
        line_index = 0
        for line in gzip.open(in2, "rt"):
            line_index += 1
            items = line.rstrip().split()
            if line_index % 4 in [1]:
                new_id = f"{items[0]}_r2rc"
                out_fh.write(f"{new_id}\n".encode('utf-8'))
                continue
            elif line_index % 4 in [2]:
                r2_seq = Seq(line.rstrip())
                r2rc_seq = r2_seq.reverse_complement()
                out_fh.write(f"{str(r2rc_seq)}\n".encode('utf-8'))
                continue
            elif line_index % 4 in [3]:
                out_fh.write("+\n".encode('utf-8'))
                continue
            elif line_index % 4 in [0]:
                out_fh.write(f"{line.rstrip()[::-1]}\n".encode('utf-8'))
                continue
            else:
                continue
        out_fh.close()

    def concat_r1_and_r2rc(self, in1):
        _cmd = list()
        _cmd.append('cat')
        _cmd.append(f"{in1}")
        _cmd.append(f"{self.r2rc_fn}")
        _cmd.append(">")
        _cmd.append(f"{self.concat_fn}")
        cmd = " ".join(_cmd)
        print(cmd)
        os.system(cmd)


def main(args):
    print(args)

    pe2se = PairedEnd2SingleEnd(args.outprefix)

    pe2se.do_r2rc(args.in2)
    pe2se.concat_r1_and_r2rc(args.in1)


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--in1", default="readyfastq/G2_1_MIA_femur_R1.fastq.gz")
    parser.add_argument("--in2", default="readyfastq/G2_1_MIA_femur_R2.fastq.gz")
    parser.add_argument("--outprefix", default="cleanfastq/G2_1_MIA_femur")
    args = parser.parse_args()
    main(args)

