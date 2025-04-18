import pandas as pd
from scipy.stats import shapiro, ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import argparse
import numpy as np
import logging

# Configure logging with time
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

def load_data(file_path):
    """Load and sort methylation data from a file."""
    data = pd.read_csv(file_path, sep='\t', low_memory=False)
    # Sort data by chrom, start, and end
    data.sort_values(by=['chrom', 'start', 'end'], inplace=True)
    return data

def process_region(region1, region2):
    """Process a single region and perform statistical tests."""
    # Extract methylation_percents
    methylation_percents1 = list(map(float, region1['methylation_percents'].split('|')))
    methylation_percents2 = list(map(float, region2['methylation_percents'].split('|')))
    
    # Check if lengths are equal
    length_mismatch = len(methylation_percents1) != len(methylation_percents2)
    
    # Initialize normality p-values
    p_value1, p_value2 = None, None
    
    # Check if data range is zero
    if np.ptp(methylation_percents1) == 0 or np.ptp(methylation_percents2) == 0:
        # Data range is zero, use Mann-Whitney U test
        test_name = "Mann-Whitney U"
        _, p_value = mannwhitneyu(methylation_percents1, methylation_percents2)
    else:
        # Perform normality test if data length is sufficient
        if len(methylation_percents1) >= 3:
            _, p_value1 = shapiro(methylation_percents1)
        if len(methylation_percents2) >= 3:
            _, p_value2 = shapiro(methylation_percents2)
        
        # Choose statistical test based on normality
        if (p_value1 is not None and p_value1 > 0.05) and (p_value2 is not None and p_value2 > 0.05):
            # Both samples are normally distributed
            test_name = "T-test"
            _, p_value = ttest_ind(methylation_percents1, methylation_percents2)
        else:
            # At least one sample is not normally distributed or normality test was not performed
            test_name = "Mann-Whitney U"
            _, p_value = mannwhitneyu(methylation_percents1, methylation_percents2)
    
    return methylation_percents1, methylation_percents2, p_value1, p_value2, test_name, p_value, length_mismatch

def analyze_files(file1, file2, sample_id_a, sample_id_b, output_file):
    """Analyze two files and write significant regions to an output file."""
    data1 = load_data(file1)
    data2 = load_data(file2)
    
    results = []
    for idx, (region1, region2) in enumerate(zip(data1.iterrows(), data2.iterrows())):
        _, region1 = region1
        _, region2 = region2
        
        # Ensure regions match
        if (region1['chrom'], region1['start'], region1['end']) != (region2['chrom'], region2['start'], region2['end']):
            logging.critical(f"Regions do not match at index {idx}: {region1['chrom']}:{region1['start']}-{region1['end']} vs {region2['chrom']}:{region2['start']}-{region2['end']}")
            raise ValueError(f"Regions do not match at index {idx}: {region1['chrom']}:{region1['start']}-{region1['end']} vs {region2['chrom']}:{region2['start']}-{region2['end']}")
        
        # Calculate distance to the next region if it exists and is on the same chromosome
        distance_to_next = None
        if idx + 1 < len(data1) and data1.iloc[idx + 1]['chrom'] == region1['chrom']:
            next_region = data1.iloc[idx + 1]
            distance_to_next = next_region['start'] - region1['end']
            distance_to_next = int(distance_to_next)

        
        methylation_percents1, methylation_percents2, p_value1, p_value2, test_name, p_value, length_mismatch = process_region(region1, region2)
        
        # Calculate mean methylation percent difference
        mean_meth_pct1 = np.mean(methylation_percents1)
        mean_meth_pct2 = np.mean(methylation_percents2)
        delta_mean_meth_pct = mean_meth_pct2 - mean_meth_pct1
        if delta_mean_meth_pct > 0:
            delta_pattern = "hyper"
        elif delta_mean_meth_pct < 0:
            delta_pattern = "hypo"
        else:
            delta_pattern = "neutral"
        
        
        # Calculate count-based methylation difference
        if length_mismatch:
            delta_cytosine_pattern = None
        elif region1['cytosine_count'] in [0]:
            delta_cytosine_pattern = None
        else:
            delta_cytosine_pattern = []
            for m1, m2 in zip(methylation_percents1, methylation_percents2):
                if m2 - m1 >= 20.0:
                    delta_cytosine_pattern.append('^')
                elif m2 - m1 <= -20.0:
                    delta_cytosine_pattern.append('_')
                else:
                    delta_cytosine_pattern.append('=')
            delta_cytosine_pattern = f"|{''.join(map(str, delta_cytosine_pattern))}|"

        results.append({
            "chrom": region1['chrom'],
            "start": region1['start'],
            "end": region1['end'],
            "distance_to_next": distance_to_next,
            "region_length": region1['region_length'],
            "cytosine_count": f"{region1['cytosine_count']}|{region2['cytosine_count']}",
            "delta_cytosine_pattern": delta_cytosine_pattern,
            f"{sample_id_a}_mean_meth_pct": round(mean_meth_pct1, 3),
            f"{sample_id_b}_mean_meth_pct": round(mean_meth_pct2, 3),
            "delta_mean_meth_pct": round(delta_mean_meth_pct, 3),
            "delta_pattern": delta_pattern,
            "significant_delta": "Y" if abs(delta_mean_meth_pct) >= 20 else "N",
            "p_value": p_value,
            "significant_p_0.01": "Y" if p_value < 0.01 else "N"
        })
    
    # FDR correction
    p_values = [result["p_value"] for result in results]
    _, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
    
    for result, corrected_p in zip(results, corrected_p_values):
        result["fdr"] = corrected_p
        result["significant_fdr_0.01"] = "Y" if corrected_p < 0.01 else "N"
    
    # Write results to output file
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_file, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze methylation data from two samples.")
    parser.add_argument("--file1", required=True, help="Path to the first methylation data file.")
    parser.add_argument("--file2", required=True, help="Path to the second methylation data file.")
    parser.add_argument("--sample_id_a", required=True, help="Sample ID for the first file.")
    parser.add_argument("--sample_id_b", required=True, help="Sample ID for the second file.")
    parser.add_argument("--output", required=True, help="Output file name for the results.")
    args = parser.parse_args()
    
    analyze_files(args.file1, args.file2, args.sample_id_a, args.sample_id_b, args.output)
