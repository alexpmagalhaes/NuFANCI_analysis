#!/usr/bin/env python

import pandas as pd
import pybedtools
from scipy.stats import poisson
import numpy as np
import argparse

def read_bedgraph(file):
    """Reads a bedgraph file and returns a pandas DataFrame."""
    return pd.read_csv(file, sep='\t', header=None, names=['chrom', 'start', 'end', 'coverage'])

def read_bed(file):
    """Reads a bed file and returns a pybedtools.BedTool object."""
    return pybedtools.BedTool(file)

def calculate_interval_coverage(bedgraph_df, interval):
    """Calculates the coverage over a specific interval from a bedgraph DataFrame."""
    chrom, start, end = interval
    subset = bedgraph_df[
        (bedgraph_df['chrom'] == chrom) &
        (bedgraph_df['end'] > start) &
        (bedgraph_df['start'] < end)
    ]
    overlap = subset.copy()
    overlap.loc[overlap['start'] < start, 'start'] = start
    overlap.loc[overlap['end'] > end, 'end'] = end

    overlap['length'] = overlap['end'] - overlap['start']
    coverage = (overlap['coverage'] * overlap['length']).sum() / (end - start)
    return coverage

def normalize_and_calculate_enrichment(input_coverage, treatment_coverage):
    """Calculates enrichment as the ratio of treatment coverage over input coverage."""
    return treatment_coverage / input_coverage if input_coverage != 0 else np.nan

def calculate_p_value(observed, expected):
    """Calculates a p-value for the observed coverage assuming a Poisson distribution with the expected mean."""
    return poisson.sf(observed, expected)

def main(input_path, test_path, bed_path):
    # Read input data
    input_bedgraph = read_bedgraph(input_path)
    treatment_bedgraph = read_bedgraph(test_path)
    intervals = read_bed(bed_path)

    # Loop through each interval in the bed file
    results = []
    counter = 0
    for interval in intervals:
        interval = interval.fields
        chrom, start, end = interval[0], int(interval[1]), int(interval[2])

        input_cov = calculate_interval_coverage(input_bedgraph, (chrom, start, end))
        treatment_cov = calculate_interval_coverage(treatment_bedgraph, (chrom, start, end))

        enrichment = normalize_and_calculate_enrichment(input_cov, treatment_cov)
        p_value = calculate_p_value(treatment_cov, input_cov)

        # Name each interval
        interval_name = f"NAD_{counter}"
        counter += 1

        results.append((chrom, start, end, interval_name, input_cov, treatment_cov, enrichment, p_value))

    # Create a DataFrame for the results
    results_df = pd.DataFrame(results, columns=['chrom', 'start', 'end', 'name', 'input_coverage', 'treatment_coverage', 'enrichment', 'p_value', 'log_p_value'])

    # Output the results with interval names
    results_df[['chrom', 'start', 'end', 'name']].to_csv('intervals_with_names.bed', index=False, sep='\t', header=False)

    # Output results for enrichment and p-values
    results_df.to_csv('coverage_enrichment_results.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate coverage enrichment over intervals with p-values.')
    parser.add_argument('--input', type=str, required=True, help='Path to the input bedgraph file.')
    parser.add_argument('--test', type=str, required=True, help='Path to the test bedgraph file.')
    parser.add_argument('--bed', type=str, required=True, help='Path to the bed file with intervals.')

    args = parser.parse_args()

    main(args.input, args.test, args.bed)
