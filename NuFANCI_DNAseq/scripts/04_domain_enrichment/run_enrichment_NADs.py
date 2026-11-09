import pandas as pd
import numpy as np
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
import argparse
import sys

def run_enrichment_analysis(input_file, output_file, coord_cols, input_cols, ip_cols):
    """
    Performs enrichment analysis on sequencing data using a Poisson model.

    This script calculates Log2 Fold Change, performs a one-sided Poisson test
    for enrichment, and corrects for multiple testing using the Benjamini-Hochberg
    method (FDR).

    Args:
        input_file (str): Path to the input CSV file.
        output_file (str): Path to save the output CSV file.
        coord_cols (list): List of column names for coordinates (e.g., ['chr', 'start', 'end']).
        input_cols (list): List of column names for Input replicates.
        ip_cols (list): List of column names for IP replicates.
    """
    try:
        # --- 1. Load Data ---
        print(f"Loading data from {input_file}...")
        df = pd.read_csv(input_file, sep = "\t")

        # Validate that all specified columns exist in the dataframe
        all_cols = coord_cols + input_cols + ip_cols
        missing_cols = [col for col in all_cols if col not in df.columns]
        if missing_cols:
            print(f"Error: The following columns were not found in the input file: {', '.join(missing_cols)}")
            sys.exit(1)

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_file}")
        sys.exit(1)

    # --- 2. Data Preparation ---
    print("Preparing data...")

    # Add a pseudocount of 1 to all count columns to avoid division by zero
    df[input_cols + ip_cols] = df[input_cols + ip_cols] + 1

    # Calculate average counts for IP and Input
    df['input_mean'] = df[input_cols].mean(axis=1)
    df['ip_mean'] = df[ip_cols].mean(axis=1)

    # --- 3. Calculate Enrichment (Fold Change) ---
    print("Calculating Fold Change (without normalization)...")

    # Calculate fold change directly on the unnormalized means.
    df['fold_change'] = df['ip_mean'] / df['input_mean']


    # --- 4. Perform Poisson Test for Significance ---
    print("Performing Poisson test (without normalization)...")

    # The null hypothesis is that the IP signal is not enriched over the Input signal.
    # The rate parameter (lambda) for the Poisson test is simply the mean of the Input counts.

    # The expected mean count in IP under the null hypothesis is the input mean
    df['expected_ip_mean'] = df['input_mean']

    # The p-value is the probability of observing a count >= the actual IP count
    # given the expected rate (lambda). This is a one-sided test for enrichment.
    # We use k-1 for the CDF because we want P(X >= k) which is 1 - P(X <= k-1).
    # We round the means because Poisson works on integer counts.
    p_values = 1 - poisson.cdf(
        k=np.round(df['ip_mean']) - 1,
        mu=df['expected_ip_mean']
    )
    df['p_value'] = p_values

    # --- 5. Multiple Testing Correction ---
    print("Adjusting p-values for multiple testing (FDR)...")

    # Use Benjamini-Hochberg method to get adjusted p-values (q-values)
    reject, q_values, _, _ = multipletests(df['p_value'], alpha=0.05, method='fdr_bh')
    df['fdr'] = q_values

    # --- 6. Save Results ---
    print(f"Saving results to {output_file}...")

    # Select columns for the final output file
    output_cols = coord_cols + ['input_mean', 'ip_mean', 'fold_change', 'p_value', 'fdr']
    final_df = df[output_cols]

    final_df.to_csv(output_file, index=False)
    print("Analysis complete!")


if __name__ == '__main__':
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description="Perform enrichment analysis on sequencing data using a Poisson model.")

    parser.add_argument('-i', '--input', required=True, help="Path to the input CSV file.")
    parser.add_argument('-o', '--output', required=True, help="Path for the output results CSV file.")
    parser.add_argument('--coords', required=True, nargs='+', help="Names of the coordinate columns (e.g., --coords chr start end).")
    parser.add_argument('--inputs', required=True, nargs='+', help="Names of the Input replicate columns (e.g., --inputs input_1 input_2 input_3).")
    parser.add_argument('--ips', required=True, nargs='+', help="Names of the IP replicate columns (e.g., --ips ip_1 ip_2 ip_3).")

    args = parser.parse_args()

    # Run the main analysis function
    run_enrichment_analysis(
        input_file=args.input,
        output_file=args.output,
        coord_cols=args.coords,
        input_cols=args.inputs,
        ip_cols=args.ips
    )
