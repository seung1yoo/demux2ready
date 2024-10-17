import hashlib
from pathlib import Path
from datetime import datetime
import argparse
import os
import logging

# 로깅 설정
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

def calculate_md5(file_path):
    """Calculate the MD5 checksum of a file."""
    logging.info(f"Calculating MD5 for: {file_path}")
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def load_existing_md5(file_path):
    """Load existing MD5 checksum from a file."""
    md5_dict = {}
    if file_path.exists():
        with open(file_path, "r") as f:
            for line in f:
                md5, path = line.strip().split()
                md5_dict[path] = md5
    return md5_dict

def main(path: Path, date: str, existing_files: list):
    # Define the file extensions to look for
    #extensions = [".cram", ".cnv.vcf.gz", ".sv.vcf.gz", ".hard-filtered.gvcf.gz", ".hard-filtered.vcf.gz", ".html"]
    extensions = [".html"]
    
    # Load existing MD5 checksums
    existing_md5 = {}
    for file in existing_files:
        existing_md5.update(load_existing_md5(Path(file)))

    # Prepare the output file
    output_file = path / f"md5sum.{date}.tsv"
    with open(output_file, "w") as out_f:
        for ext in extensions:
            for file in path.rglob(f"*{ext}"):
                file_path_str = str(file)
                if file_path_str in existing_md5:
                    md5 = existing_md5[file_path_str]
                else:
                    md5 = calculate_md5(file)
                out_f.write(f"{md5}\t{file_path_str}\n")
                print(f"Processed: {file_path_str}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=Path, default=Path.cwd())
    parser.add_argument("--date", type=str, default=datetime.now().strftime("%Y-%m-%d"))
    parser.add_argument("--existing", nargs="+", help="List of files with existing MD5 checksums", default=[])
    args = parser.parse_args()
    
    # Convert args.path to an absolute path
    args.path = args.path.resolve()
    
    main(args.path, args.date, args.existing)
