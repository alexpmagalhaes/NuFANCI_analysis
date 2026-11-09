"""Step 00: filter the FragPipe protein table before differential analysis.

Removes from FragPipe's combined_protein.tsv:
  1. proteins whose "Protein ID" is in the contaminant list (reference_lists/Contaminants.tsv);
  2. proteins named "contam_sp..." or containing "Cont_" (contaminant entries in the search FASTA);
  3. proteins with "Combined Total Peptides" <= 1 (single-peptide identifications).
Writes <input>_f.tsv (e.g. combined_protein_f.tsv), the input to FragPipe-Analyst and the notebooks.

Usage:
    python filter_proteins.py -i combined_protein.tsv -c ../reference_lists/Contaminants.tsv
"""
import pandas as pd
import argparse
import os

def filter_contaminants(input_file, contaminants_file):
    """
    Filters rows from a master TSV file based on a list of Protein IDs,
    protein names starting with 'contam_sp' or containing 'sp|Cont_',
    and the number of combined total peptides.

    Args:
        input_file (str): The path to the master data file.
        contaminants_file (str): The path to the file with contaminant Protein IDs.
    """
    try:
        # --- 1. Load the contaminant proteins ---
        # Read the TSV file containing the contaminant protein IDs.
        # We expect a column named "Protein_ID" in this file.
        print(f"Reading contaminants from: {contaminants_file}")
        contaminants_df = pd.read_csv(contaminants_file, sep='\t')

        # Check if the required column exists
        if 'Protein_ID' not in contaminants_df.columns:
            print(f"Error: The contaminants file '{contaminants_file}' must contain a column named 'Protein_ID'.")
            return

        # Create a set of contaminant protein IDs for fast lookups.
        contaminant_ids = set(contaminants_df['Protein_ID'])
        print(f"Found {len(contaminant_ids)} unique contaminant protein IDs.")

        # --- 2. Load the master data file ---
        # Read the main input TSV file.
        # We expect columns named "Protein ID", "Protein", and "Combined Total Peptides".
        print(f"Reading master data from: {input_file}")
        master_df = pd.read_csv(input_file, sep='\t')

        # Check if the required columns exist
        if 'Protein ID' not in master_df.columns:
            print(f"Error: The input file '{input_file}' must contain a column named 'Protein ID'.")
            return
        if 'Protein' not in master_df.columns:
            print(f"Error: The input file '{input_file}' must contain a column named 'Protein'.")
            return
        if 'Combined Total Peptides' not in master_df.columns:
            print(f"Error: The input file '{input_file}' must contain a column named 'Combined Total Peptides'.")
            return

        initial_rows = len(master_df)
        print(f"Master table contains {initial_rows} rows.")

        # --- 3. Filter the data ---
        # Condition 1: Keep rows where "Protein ID" is NOT in the contaminant_ids set.
        id_filter = ~master_df['Protein ID'].isin(contaminant_ids)

        # Condition 2: Keep rows where "Protein" name does NOT start with "contam_sp"
        # and does NOT contain "sp|Cont_" (another naming convention used for contaminants).
        # The .str accessor works on Series with string values. na=False treats NaN as not matching.
        name_filter = (
            ~master_df['Protein'].str.startswith('contam_sp', na=False)
            & ~master_df['Protein'].str.contains('Cont_', regex=False, na=False)
        )

        # Condition 3: Keep rows where "Combined Total Peptides" is greater than 1.
        # First, ensure the column is treated as a numeric type. `errors='coerce'` will turn non-numeric values into NaN.
        # Then, we fill any potential NaN values with 0 so the comparison works correctly.
        peptides_column = pd.to_numeric(master_df['Combined Total Peptides'], errors='coerce').fillna(0)
        peptide_filter = peptides_column > 1

        # Combine the filters: a row is kept only if ALL three conditions are true.
        print("Filtering rows based on contaminant list, 'contam_sp'/'sp|Cont_' naming, and peptide count...")
        filtered_df = master_df[id_filter & name_filter & peptide_filter]

        final_rows = len(filtered_df)
        rows_removed = initial_rows - final_rows
        print(f"Removed {rows_removed} rows. The filtered table has {final_rows} rows.")

        # --- 4. Save the output file ---
        # Construct the output filename by adding "_f" before the extension.
        base, ext = os.path.splitext(input_file)
        output_file = f"{base}_f{ext}"

        # Save the filtered DataFrame to a new TSV file without the pandas index.
        filtered_df.to_csv(output_file, sep='\t', index=False)

        print(f"\nSuccessfully created filtered file: {output_file}")

    except FileNotFoundError as e:
        print(f"Error: File not found - {e.filename}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    # Set up the command-line argument parser
    parser = argparse.ArgumentParser(
        description="Filter a master TSV file by removing contaminants and rows with low peptide counts."
    )

    # Define command-line arguments
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to the input master TSV file.'
    )
    parser.add_argument(
        '-c', '--contaminants',
        required=True,
        help='Path to the contaminants TSV file (must contain a "Protein_ID" column).'
    )

    # Parse the arguments provided by the user
    args = parser.parse_args()

    # Run the main filtering function
    filter_contaminants(args.input, args.contaminants)
