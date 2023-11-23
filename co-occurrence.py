#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np

# Argument parser
parser = argparse.ArgumentParser(description='Correlation analysis script')
parser.add_argument('-i', '--input', type=str, required=True,
                    help='Abundance table file (tsv)')
parser.add_argument('-o', '--output', type=str, required=True,
                    help='Prefix for the output name.')
parser.add_argument('-p', '--prevalence', type=float, default=0.1,
                    help='Minimum percentage of samples for correlation analysis (default: 0.1)')
parser.add_argument('-c', '--correlation', type=float, default=0.3,
                    help='Minimum correlation to keep pairs (default: 0.3)')
args = parser.parse_args()

print("Read in abundance table.")
OTU = pd.read_csv(args.input, sep="\t", index_col=0)

df = OTU.drop(OTU.columns[OTU.columns.str.contains('NC')], axis=1)
df['sample_count'] = df.apply(lambda row: row[row != 0].count(), axis=1)
df['proportion_samples'] = df['sample_count'] / (df.shape[1] - 1)

# Define the threshold
threshold = args.prevalence

# Filter rows where the proportion of 0s is less than or equal to the threshold
filtered_df = df[df['proportion_samples'] >= threshold].drop(['sample_count', 'proportion_samples'], axis=1)

prevalence_df = filtered_df.map(lambda x: 1 if x > 0 else x)

n = len(prevalence_df)

print(f"Calculate correlation matrix for {n} contigs (contig prevalence in samples = {threshold*100}%).")
prevalence_df_transposed = prevalence_df.transpose()
correlation_matrix=prevalence_df_transposed.corr()

mask = np.triu(np.ones(correlation_matrix.shape), k=1).astype(bool)
masked_correlation_matrix = correlation_matrix.mask(mask)
correlation_matrix=masked_correlation_matrix.rename_axis(axis=0, mapper="Contig1").rename_axis(axis=1, mapper="Contig2")

print("Write correlation matrix.")
correlation_matrix.to_csv(args.output+"_correlation_matrix.tsv", sep='\t', index=True)

print("Write pairwise dataframe.")

related_contigs=correlation_matrix[abs(correlation_matrix) >= args.correlation].stack()

result_df = pd.DataFrame(related_contigs)
result_df = result_df.reset_index()

# Rename existing columns if necessary
result_df.columns = ["Contig1", "Contig2", "Correlation"]
result_df = result_df[result_df['Contig1'] != result_df['Contig2']]

result_df.sort_values(by='Contig1', inplace=True)
result_df.to_csv(args.output+"_related_contigs.tsv", sep='\t', index=False)

print("Finished.")

