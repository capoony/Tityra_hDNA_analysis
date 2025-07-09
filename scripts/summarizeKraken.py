import pandas as pd
import argparse
import glob
import os

# Function to parse Kraken reports with the described whitespace handling and hierarchical logic


def parse_kraken_report_hierarchy(file_path):
    data = []
    last_kingdom = None  # To track the last Kingdom-level taxon
    last_order = None    # To track the last Order-level taxon

    with open(file_path, 'r') as f:
        for line in f:
            # Remove leading spaces and normalize tabs
            line = line.lstrip()
            parts = line.split('\t')  # Split by tabs

            # Check if the line has enough fields to process
            if len(parts) >= 6:
                try:
                    # % abundance in the first column
                    abundance = float(parts[0])
                except ValueError:
                    continue  # Skip lines with invalid numeric values

                taxon_name = parts[-1].strip()  # Taxon name in the last column
                # Hierarchical level (e.g., K, O, G)
                taxon_level = parts[3].strip()

                # Retain the last Kingdom and Order names
                if taxon_level == "D":
                    last_kingdom = taxon_name
                elif taxon_level == "O":
                    last_order = taxon_name

                # Process Genus-level taxa
                if taxon_level == "S" and abundance >= 0.05:
                    if last_kingdom and last_order:  # Ensure Kingdom and Order are not None
                        full_taxon = f"{last_kingdom}; {last_order}; {taxon_name}"
                        data.append([full_taxon, abundance])

    return pd.DataFrame(data, columns=["Taxon", "Abundance (%)"])

# Main function to handle input and output using argparse


def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(
        description="Process Kraken reports and generate a merged table.")
    parser.add_argument("--input", required=True,
                        help="Input folder containing Kraken report files.")
    parser.add_argument("--output", required=True,
                        help="Output CSV file to save the merged table.")
    args = parser.parse_args()

    # Get all Kraken report files from the input directory
    report_files = glob.glob(f"{args.input}/*.txt")

    if not report_files:
        print(f"No report files found in {args.input}.")
        return

    # Parse all Kraken report files and extract sample names from filenames
    parsed_reports = {}
    for file_path in report_files:
        sample_name = os.path.basename(file_path).rsplit(".txt", 1)[0]
        parsed_reports[sample_name] = parse_kraken_report_hierarchy(file_path)

    # Merge all parsed reports into a single table, filling missing taxa with NA
    merged_data = pd.DataFrame()
    for sample_name, df in parsed_reports.items():
        if merged_data.empty:
            merged_data = df.set_index("Taxon")
            merged_data.rename(
                columns={"Abundance (%)": sample_name}, inplace=True)
        else:
            merged_data = merged_data.join(
                df.set_index("Taxon").rename(
                    columns={"Abundance (%)": sample_name}),
                how="outer"
            )

    # Fill missing values with NA
    merged_data = merged_data.reset_index()

    # Save the merged data to the specified output file
    merged_data.to_csv(args.output, index=False)
    print(f"Merged Kraken data saved to {args.output}.")


if __name__ == "__main__":
    main()
