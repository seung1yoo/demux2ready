import hashlib
from pathlib import Path
from datetime import datetime
import argparse
import os
import logging
from concurrent.futures import ProcessPoolExecutor

# 로깅 설정
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

def calculate_md5_and_size(file_path):
    """Calculate the MD5 checksum and size of a file."""
    logging.info(f"Calculating MD5 and size for: {file_path}")
    hash_md5 = hashlib.md5()
    total_size = 0
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
            total_size += len(chunk)
    return hash_md5.hexdigest(), total_size, str(file_path.resolve())

def process_file(file, existing_md5):
    file_path_str = str(file.resolve())
    if file_path_str in existing_md5:
        md5 = existing_md5[file_path_str][0]
        size = existing_md5[file_path_str][1]
    else:
        md5, size, file_path_str = calculate_md5_and_size(file)
    return md5, size, file_path_str

def load_existing_md5(file_path):
    """Load existing MD5 checksums from a file."""
    existing_md5 = {}
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")

            if len(parts) in [3]:
                md5, size, path = parts
            elif len(parts) in [2]:
                logging.warning(f"Invalid format in line: {line.strip()} -> the file size is not provided. calculating...")
                md5, path = parts
                size = os.path.getsize(path)
            else:
                logging.warning(f"Invalid format in line: {line.strip()} -> skipping")
                continue

            try:
                size = int(size)
            except ValueError:
                logging.warning(f"Invalid size in line: {line.strip()}")
                continue
            existing_md5[path] = (md5, size)
    return existing_md5

def main(path: Path, date: str, existing_files: list, num_cpus: int):
    # Define the file extensions to look for
    extensions = [".cram", ".cnv.vcf.gz", ".sv.vcf.gz", ".hard-filtered.gvcf.gz", ".hard-filtered.vcf.gz", ".html"]
    
    # Load existing MD5 checksums
    existing_md5 = {}
    for file in existing_files:
        existing_md5.update(load_existing_md5(Path(file)))

    # Prepare the output file
    output_file = path / f"md5sum.{date}.tsv"
    with open(output_file, "w") as out_f:
        with ProcessPoolExecutor(max_workers=num_cpus) as executor:
            futures = [executor.submit(process_file, file, existing_md5) for ext in extensions for file in path.rglob(f"*{ext}")]
            for future in futures:
                md5, size, file_path_str = future.result()
                out_f.write(f"{md5}\t{size}\t{file_path_str}\n")
                print(f"Processed: {file_path_str}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=Path, default=Path.cwd())
    parser.add_argument("--date", type=str, default=datetime.now().strftime("%Y-%m-%d"))
    parser.add_argument("--existing", nargs="+", help="List of files with existing MD5 checksums", default=[])
    parser.add_argument("--cpus", type=int, default=1, help="Number of CPUs to use for parallel processing")
    args = parser.parse_args()
    
    # Convert args.path to an absolute path
    args.path = args.path.resolve()
    
    main(args.path, args.date, args.existing, args.cpus)
