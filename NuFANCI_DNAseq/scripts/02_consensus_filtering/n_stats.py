import sys

def calculate_n_stats(file_path, n_percentages=[5, 10, 25, 50, 75, 90, 95]):
    """
    Calculates N-statistics for intervals in a BED file.

    Args:
        file_path (str): The path to the BED file.
        n_percentages (list): A list of percentages for which to calculate N-statistics.

    Returns:
        dict: A dictionary containing the total size and the calculated N-statistics.
    """

    try:
        with open(file_path, 'r') as f:
            lengths = []
            for line in f:
                # Skip header lines if they exist
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')

                # A valid BED line has at least 3 columns
                if len(parts) >= 3:
                    try:
                        start = int(parts[1])
                        end = int(parts[2])
                        length = end - start

                        if length > 0:
                            lengths.append(length)
                    except (ValueError, IndexError):
                        # Skip malformed lines
                        continue

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        sys.exit(1)

    if not lengths:
        return {"Total Size": 0, "N-stats": {f"N{p}": 0 for p in n_percentages}}

    # Sort the lengths in descending order
    lengths.sort(reverse=True)

    total_len = sum(lengths)
    n_stats = {}

    # Calculate thresholds for each N-statistic
    thresholds = {p: total_len * (p / 100.0) for p in n_percentages}

    current_sum = 0

    for length in lengths:
        current_sum += length

        # Check if the current sum has met or exceeded any threshold
        for p in n_percentages:
            if f"N{p}" not in n_stats and current_sum >= thresholds[p]:
                n_stats[f"N{p}"] = length

    # Handle cases where some N-stats might not be reached (e.g., if total length is very small)
    for p in n_percentages:
        if f"N{p}" not in n_stats:
            n_stats[f"N{p}"] = 0

    return {"Total Size": total_len, "N-stats": n_stats}

# --- Main execution block ---
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <bed_file_path>")
        sys.exit(1)

    bed_file = sys.argv[1]
    results = calculate_n_stats(bed_file)

    print(f"Total Size: {results['Total Size']}")
    for stat, value in results['N-stats'].items():
        print(f"{stat}: {value}")
