#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# Argument parser
parser = argparse.ArgumentParser(
    description="Correlation analysis script for specific contigs."
)
parser.add_argument(
    "-i", "--input", type=str, required=True, help="Abundance table file (tsv)"
)
parser.add_argument(
    "-s",
    "--segments",
    type=str,
    required=True,
    help="File with a list of contigs of interest (often RdRP segments), each on a new line.",
)
parser.add_argument(
    "-o", "--output", type=str, required=True, help="Prefix for the output name."
)
parser.add_argument(
    "-p",
    "--prevalence",
    type=float,
    default=0.1,
    help="Minimum percentage of samples for correlation analysis (default: 0.1)",
)
parser.add_argument(
    "-c",
    "--correlation",
    type=float,
    default=0.3,
    help="Minimum correlation to keep pairs (default: 0.3)",
)
parser.add_argument(
    "-l",
    "--lengths",
    type=str,
    help="File with the lengths of each contig",
)
args = parser.parse_args()

print("Read in abundance table.")
OTU = pd.read_csv(args.input, sep="\t", index_col=0)

df = OTU.drop(OTU.columns[OTU.columns.str.contains("NC")], axis=1)
df["sample_count"] = df.apply(lambda row: row[row != 0].count(), axis=1)
df["proportion_samples"] = df["sample_count"] / (df.shape[1] - 1)

# Define the threshold
threshold = args.prevalence

# Filter rows where the proportion of 0s is less than or equal to the threshold
filtered_df = df[df["proportion_samples"] >= threshold].drop(
    ["sample_count", "proportion_samples"], axis=1
)

n = len(filtered_df)

print(
    f"{n} contigs were retained for correlation analysis (prevalence in samples = {threshold*100}%)"
)

if args.lengths:
    lengths = pd.read_csv(
        args.lengths, sep="\t", index_col=0, header=None, names=["Contig", "length"]
    )
    df = filtered_df.div(lengths["length"], axis=0).dropna(how="all")
else:
    df = filtered_df.map(lambda x: 1 if x > 0 else x)

file_path = Path(args.segments)

# Read the file into a list
with open(file_path, "r") as file:
    segment_list = file.readlines()

segment_list = [line.strip() for line in segment_list]

# Create an empty DataFrame to store the results
correlation_results_df = pd.DataFrame()

# Loop through each segment in segment_list
for i in segment_list:
    df2 = df.loc[i]
    df3 = df.corrwith(df2, axis=1, method="spearman")

    # Append the results to the DataFrame with i as the column name
    correlation_results_df[i] = df3

cor_threshold = args.correlation

mask = (abs(correlation_results_df) >= cor_threshold).any(axis=1)
corr_df = correlation_results_df[mask]

print(f"Write correlation matrix with a threshold of {cor_threshold}")
corr_df.to_csv(args.output + ".tsv", sep="\t", index=True)

print("Finished.")
