# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 14:25:09 2024

@author: pahermans
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV file into a DataFrame
algo = "SVC"
df = pd.read_csv("E:/VM_data/Git/VCF_compression_benchmarking/Data/"+algo+"_ALLchrom.csv")
#df = pd.read_csv("E:/VM_data/Git/VCF_compression_benchmarking/Data/SVC_ALLchrom.csv")

# Filter the rows by Command to get separate DataFrames for 'compress' and 'decompress'
compress_df = df[df['Command'] == 'compress'][['Chromosome', 'Sample_Size', 'Input_file_size']].groupby(['Chromosome', 'Sample_Size'], as_index=False).mean()
decompress_df = df[df['Command'] == 'decompress'][['Chromosome', 'Sample_Size', 'Input_file_size']].groupby(['Chromosome', 'Sample_Size'], as_index=False).mean()



# Merge the DataFrames on Chromosome and Sample_Size to align compress and decompress values
merged_df = pd.merge(compress_df, decompress_df, on=['Chromosome', 'Sample_Size'], suffixes=('_compress', '_decompress'))

# Calculate the Compression_Ratio
merged_df['Compression_Ratio'] = merged_df['Input_file_size_compress'] / merged_df['Input_file_size_decompress']

# Normalize the Compression_Ratio by the Input_file_size of the compress command
#merged_df['Normalized_Compression_Ratio'] = merged_df['Compression_Ratio'] / merged_df['Input_file_size_compress']

# Select only the necessary columns
compression_ratio_df = merged_df[['Chromosome', 'Sample_Size', 'Compression_Ratio', 'Input_file_size_compress']]
# From KB to GB
compression_ratio_df.loc[:,'Input_file_size_compress'] = compression_ratio_df['Input_file_size_compress'] / (1024**3)

# Define the colormap with unique colors
num_colors = 22
cmap = plt.cm.get_cmap("gist_rainbow", num_colors)  # 'tab20' or 'hsv', or any other suitable colormap
# Assign colors to each chromosome
color_map = {chromosome: cmap(chromosome) for i, chromosome in enumerate(np.sort(df['Chromosome'].unique()))}

# Display the result DataFrame
print(compression_ratio_df)


# Plot Compression_Ratio vs Sample_Size with a line for each chromosome
plt.figure(figsize=(10, 6))

# Group by Chromosome to plot each one separately
for chromosome, group in compression_ratio_df.groupby('Chromosome'):
    plt.plot(group['Input_file_size_compress'], group['Compression_Ratio'], marker='o', color=color_map[chromosome], label=f'Chr {chromosome}')


# Add labels and title
plt.xlabel('Input file size (uncompressed) in GB')
plt.ylabel('Compression Ratio')
#plt.xscale("log")
titre = algo +' - Compression Ratio vs Input filesize for Each Chromosome'
plt.title(titre)
plt.legend(title='Chromosomes',bbox_to_anchor=(1.2, 1.1))
plt.grid(True)

plt.savefig("E:\VM_data\Git\VCF_compression_benchmarking\Images\\"+titre , bbox_inches='tight')
# Display the plot
plt.show()


plt.figure(figsize=(10, 6))
# Group by Chromosome to plot each one separately
for chromosome, group in compression_ratio_df.groupby('Chromosome'):
    plt.plot(group['Sample_Size'], group['Compression_Ratio'], marker='o', color=color_map[chromosome], label=f'Chr {chromosome}')


# Add labels and title
plt.xlabel('Sample Size')
plt.ylabel('Compression Ratio')
#plt.xscale("log")
titre = algo +' - Compression Ratio vs Sample Size for Each Chromosome'
plt.title(titre)
plt.legend(title='Chromosomes',bbox_to_anchor=(1.2, 1.1))
plt.grid(True)

plt.savefig("E:\VM_data\Git\VCF_compression_benchmarking\Images\\"+titre ,bbox_inches='tight' )

# Display the plot
plt.show()