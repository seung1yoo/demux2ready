import logging
import pandas as pd
import gzip
import argparse
from pathlib import Path
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s'
)


def load_methylation_data(chromosome, sample, context, home_dir):
    filename = f"{sample}_R1_bismark_bt2_pe.deduplicated.{context}_report.{chromosome}.txt.gz"
    filepath = Path(home_dir) / "bismark" / f"{sample}_split_cx_to_context" / filename

    logging.info(f"Loading methylation data from {filepath}")

    if not filepath.exists():
        logging.warning(f"File {filepath} not found.")
        return pd.DataFrame()
    
    with gzip.open(filepath, 'rt') as f:
        methylation_data = pd.read_csv(f, sep='\t', header=None, 
                                       names=['chrom', 'position', 'strand', 'meth_count', 'unmeth_count', 'context', 'sequence'])
    return methylation_data

def process_chromosome(chrom, group, sample, context, min_depth, min_cytosine, home_dir):
    logging.info(f"Processing chromosome {chrom}")
    methylation_data = load_methylation_data(chrom, sample, context, home_dir)
    logging.info(f"Finished loading methylation data for chromosome {chrom}")
    if methylation_data.empty:
        logging.info(f"No methylation data found for chromosome {chrom}")
        return []

    chrom_results = []
    for _, region in group.iterrows():
        region_data = methylation_data[(methylation_data['position'] >= region['start']) & 
                                       (methylation_data['position'] <= region['end'])].copy()
        region_data.loc[:, 'depth'] = region_data['meth_count'] + region_data['unmeth_count']
        region_data = region_data[region_data['depth'] >= min_depth]

        region_length = region['end'] - region['start'] + 1

        if not region_data.empty and len(region_data) >= min_cytosine:
            methylation_percents = np.where(
                region_data['depth'] == 0, 0.0, 
                region_data['meth_count'] / region_data['depth'] * 100
            ).round(2).astype(str).tolist()
            meth_counts = region_data['meth_count'].astype(str).tolist()
            unmeth_counts = region_data['unmeth_count'].astype(str).tolist()
            depth_s = region_data['depth'].astype(str).tolist()
            cytosine_count = len(region_data)
        else:
            methylation_percents = ['0']
            meth_counts = ['0']
            unmeth_counts = ['0']
            depth_s = ['0']
            cytosine_count = 0

        chrom_results.append([
            region['chrom'], region['start'], region['end'], region_length, 
            cytosine_count, '|'.join(methylation_percents), 
            '|'.join(meth_counts), '|'.join(unmeth_counts), '|'.join(depth_s)
        ])


    num_target_regions = len(chrom_results)
    avg_region_length = np.mean([result[3] for result in chrom_results])
    avg_cytosine_count = np.mean([result[4] for result in chrom_results if result[4] > min_cytosine])
    depth_values = [int(d) for result in chrom_results for d in result[8].split('|') if int(d) >= min_depth]
    avg_depth = np.mean(depth_values)

    logging.info(f"Chromosome {chrom} statistics: "
                 f"Number of target regions: {num_target_regions}, "
                 f"Average region length: {avg_region_length:.2f}, "
                 f"Average cytosine count: {avg_cytosine_count:.2f}, "
                 f"Average depth: {avg_depth:.2f}")

    return chrom_results

def process_target_regions(target_file, sample, context, min_depth, min_cytosine, home_dir):
    logging.info(f"Processing target regions from {target_file}")
    target_regions = pd.read_csv(target_file, sep='\t', header=None, names=['chrom', 'start', 'end'], low_memory=False)
    results = []

    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_chromosome, chrom, group, sample, context, min_depth, min_cytosine, home_dir): chrom 
                   for chrom, group in target_regions.groupby('chrom')}
        
        for future in as_completed(futures):
            chrom_results = future.result()
            results.extend(chrom_results)

    results_df = pd.DataFrame(results, columns=[
        'chrom', 'start', 'end', 'region_length', 'cytosine_count', 
        'methylation_percents', 'meth_counts', 'unmeth_counts', 'depth_s'
    ])

    results_df.sort_values(by=['chrom', 'start'], inplace=True)

    output_file = (
        Path(home_dir) / "bismark" / f"{sample}_split_cx_to_context" / 
        f"{sample}_R1_bismark_bt2_pe.deduplicated.targetRegion.{context}.minD{min_depth}.minC{min_cytosine}.bed"
    )
    logging.info(f"Writing sorted results to {output_file}")
    results_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract methylation information for target regions.")
    parser.add_argument('--target', required=True, help='Path to the Target_region.bed file')
    parser.add_argument('--home_dir', required=True, help='Directory where methylation data files are located')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--context', required=True, help='Cytosine context (e.g., CG, CHG, CHH)')
    parser.add_argument('--min_depth', type=int, required=True, help='Minimum depth on A cytosine base for filtering')
    parser.add_argument('--min_cytosine', type=int, required=True, help='Minimum number of cytosines on A target region for filtering')
    args = parser.parse_args()

    logging_filename = f"bismark_cx_to_targetRegion_{args.sample}_{args.context}_minD{args.min_depth}_minC{args.min_cytosine}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    console = logging.FileHandler(filename=logging_filename, mode='w')
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

    logging.info(f"Arguments: {args}")
    process_target_regions(args.target, args.sample, args.context, args.min_depth, args.min_cytosine, args.home_dir)



