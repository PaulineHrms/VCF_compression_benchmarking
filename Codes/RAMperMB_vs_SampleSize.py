# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 17:01:12 2024

@author: pahermans
"""

import pandas as pd
import matplotlib.pyplot as plt


#%%
################# ALL CHROM, SOME TESTS ##########################
# Load the CSV file
algo = "GSC"
df = pd.read_csv("E:/VM_data/Git/VCF_compression_benchmarking/Data/"+algo+"_ALLchrom.csv")


# Convert Input_file_size from bytes to MB
df['Input_file_size_MB'] = df['Input_file_size'] / 1048576

# Calculate Time per MB
df['RAM_per_MB'] = df['Memory_Used_MB'] / df['Input_file_size_MB']

# Filter data for each command
compress_df = df[df['Command'] == 'compress'][['Chromosome', 'Sample_Size',  'Input_file_size_MB', 'RAM_per_MB','Memory_Used_MB']].groupby(['Chromosome', 'Sample_Size'], as_index=False).mean()
decompress_df = df[df['Command'] == 'decompress'][['Chromosome', 'Sample_Size',  'Input_file_size_MB', 'RAM_per_MB','Memory_Used_MB']].groupby(['Chromosome', 'Sample_Size'], as_index=False).mean()

# For random_access, calculate the mean Time_per_MB for each Sample_Size
random_access_df = df[df['Command'] == 'random_access'][['Chromosome', 'Sample_Size', 'Input_file_size_MB', 'RAM_per_MB','Memory_Used_MB']].groupby(['Chromosome', 'Sample_Size'], as_index=False).mean()

ALL_df = df[['Command','Sample_Size',  'Input_file_size_MB', 'RAM_per_MB','Memory_Used_MB']].groupby(['Command','Sample_Size'], as_index=False).mean()
ALL_compress_df = ALL_df[ALL_df['Command'] == 'compress']
ALL_decompress_df = ALL_df[ALL_df['Command'] == 'decompress']
ALL_random_access_df = ALL_df[ALL_df['Command'] == 'random_access']

unique_chromosomes = df['Chromosome'].unique()
# Define the colormap with unique colors
num_colors = 22
cmap = plt.cm.get_cmap("gist_rainbow", num_colors)  # 'tab20' or 'hsv', or any other suitable colormap
# Assign colors to each chromosome
color_map = {chromosome: cmap(chromosome) for i, chromosome in enumerate(unique_chromosomes)}


# Plot for the "compress" command
plt.figure(figsize=(8, 5))
for chromosome, group in compress_df.groupby('Chromosome'):
    plt.plot(group['Sample_Size'], group['RAM_per_MB'], marker='o', color=color_map[chromosome],  label=f'Chr {chromosome}')
plt.plot(ALL_compress_df['Sample_Size'], ALL_compress_df['RAM_per_MB'],marker='o', color='black', label="Mean")
plt.xlabel('Sample Size')
plt.ylabel('RAM per MB ')
plt.xscale("log")
titre = algo +' - RAM per MB vs Sample Size for Compress Command'
plt.title(titre)
plt.legend(title='Chromosomes',bbox_to_anchor=(1.2, 1.1))
plt.grid(True)
plt.savefig("E:\VM_data\Git\VCF_compression_benchmarking\Images\\" + titre, bbox_inches='tight')

plt.show()

# Plot for the "decompress" command
plt.figure(figsize=(8, 5))
for chromosome, group in decompress_df.groupby('Chromosome'):
    plt.plot(group['Sample_Size'], group['RAM_per_MB'], marker='o', color=color_map[chromosome], label=f'Chr {chromosome}')
plt.plot(ALL_decompress_df['Sample_Size'], ALL_decompress_df['RAM_per_MB'],marker='o', color='black', label="Mean")
plt.xlabel('Sample Size')
plt.ylabel('RAM per MB')
plt.xscale("log")
titre = algo +' - RAM per MB vs Sample Size for Decompress Command'
plt.title(titre)
plt.legend(title='Chromosomes',bbox_to_anchor=(1.2, 1.1))
plt.grid(True)
plt.savefig("E:\VM_data\Git\VCF_compression_benchmarking\Images\\" + titre, bbox_inches='tight')
plt.show()

# Plot for the "random_access" command with averaged Time_per_MB
plt.figure(figsize=(8, 5))
for chromosome, group in random_access_df.groupby('Chromosome'):
    plt.plot(group['Sample_Size'], group['RAM_per_MB'], marker='o', color=color_map[chromosome],  label=f'Chr {chromosome}')
plt.plot(ALL_random_access_df['Sample_Size'], ALL_random_access_df['RAM_per_MB'],marker='o', color='black', label="Mean")
plt.xlabel('Sample Size')
plt.xscale("log")
plt.ylabel('RAM per MB')
titre = algo +' - RAM per MB vs Sample Size for Random Access Command'
plt.title(titre)
plt.legend(title='Chromosomes',bbox_to_anchor=(1.2, 1.1))
plt.grid(True)
plt.savefig("E:\VM_data\Git\VCF_compression_benchmarking\Images\\" + titre, bbox_inches='tight')
plt.show()


#%%
################# SOME CHROM, FULL TESTS ##########################

if algo=="GSC":
    # Load the CSV file
    df = pd.read_csv("E:/VM_data/Git/VCF_compression_benchmarking/Data/GSC_additional_tests.csv")
    
    
    # Convert Input_file_size from bytes to MB
    df['Input_file_size_MB'] = df['Input_file_size'] / 1048576
    
    # Calculate Time per MB
    df['RAM_per_MB'] = df['Memory_Used_MB'] / df['Input_file_size_MB']
    
    # Filter data for each command
    compress_df = df[df['Command'] == 'compress'][['Chromosome', 'Sample_Size',  'Input_file_size_MB', 'RAM_per_MB','Memory_Used_MB']].groupby(['Chromosome', 'Sample_Size'], as_index=False).mean()
    decompress_df = df[df['Command'] == 'decompress'][['Chromosome', 'Sample_Size',  'Input_file_size_MB', 'RAM_per_MB','Memory_Used_MB']].groupby(['Chromosome', 'Sample_Size'], as_index=False).mean()
    decompress_one_sample_df = df[df['Command'] == 'decompress_one_sample'][['Chromosome', 'Sample_Size', 'Input_file_size_MB', 'RAM_per_MB','Memory_Used_MB']].groupby(['Chromosome', 'Sample_Size'], as_index=False).mean()
    
    # For random_access, calculate the mean Time_per_MB for each Sample_Size
    random_access_df = df[df['Command'] == 'random_access'][['Chromosome', 'Sample_Size', 'Input_file_size_MB', 'RAM_per_MB','Memory_Used_MB']].groupby(['Chromosome', 'Sample_Size'], as_index=False).mean()
    random_access_one_sample_df = df[df['Command'] == 'random_access_one_sample'][['Chromosome', 'Sample_Size', 'Input_file_size_MB', 'RAM_per_MB','Memory_Used_MB']].groupby(['Chromosome', 'Sample_Size'], as_index=False).mean()
    
    ALL_df = df[['Command','Sample_Size',  'Input_file_size_MB', 'RAM_per_MB','Memory_Used_MB']].groupby(['Command','Sample_Size'], as_index=False).mean()
    ALL_compress_df = ALL_df[ALL_df['Command'] == 'compress']
    ALL_decompress_df = ALL_df[ALL_df['Command'] == 'decompress']
    ALL_decompress_one_sample_df = ALL_df[ALL_df['Command'] == 'decompress_one_sample']
    ALL_random_access_df = ALL_df[ALL_df['Command'] == 'random_access']
    ALL_random_access_one_sample_df = ALL_df[ALL_df['Command'] == 'random_access_one_sample']
    
    
    unique_chromosomes = df['Chromosome'].unique()
    # Define the colormap with unique colors
    num_colors = 22
    cmap = plt.cm.get_cmap("gist_rainbow", num_colors)  # 'tab20' or 'hsv', or any other suitable colormap
    # Assign colors to each chromosome
    color_map = {chromosome: cmap(chromosome) for i, chromosome in enumerate(unique_chromosomes)}
    
    
    # Plot for the "compress" command
    plt.figure(figsize=(8, 5))
    for chromosome, group in compress_df.groupby('Chromosome'):
        plt.plot(group['Sample_Size'], group['RAM_per_MB'], marker='o', color=color_map[chromosome],  label=f'Chr {chromosome}')
    plt.plot(ALL_compress_df['Sample_Size'], ALL_compress_df['RAM_per_MB'],marker='o', color='black', label="Mean")
    plt.xlabel('Sample Size')
    plt.ylabel('RAM per MB ')
    plt.xscale("log")
    titre = algo +' - RAM per MB vs Sample Size for Compress Command (bis)'
    plt.title(titre)
    plt.legend(title='Chromosomes')#,bbox_to_anchor=(1.2, 1.1))
    plt.grid(True)
    plt.savefig("E:\VM_data\Git\VCF_compression_benchmarking\Images\\" + titre, bbox_inches='tight')
    plt.show()
    
    # Plot for the "decompress" command
    plt.figure(figsize=(8, 5))
    for chromosome, group in decompress_df.groupby('Chromosome'):
        plt.plot(group['Sample_Size'], group['RAM_per_MB'], marker='o', color=color_map[chromosome], label=f'Chr {chromosome}')
    plt.plot(ALL_decompress_df['Sample_Size'], ALL_decompress_df['RAM_per_MB'],marker='o', color='black', label="Mean")
    plt.xlabel('Sample Size')
    plt.ylabel('RAM per MB')
    plt.xscale("log")
    titre = algo +' - RAM per MB vs Sample Size for Decompress Command (bis)'
    plt.title(titre)
    plt.legend(title='Chromosomes')#,bbox_to_anchor=(1.2, 1.1))
    plt.grid(True)
    plt.savefig("E:\VM_data\Git\VCF_compression_benchmarking\Images\\" + titre, bbox_inches='tight')
    plt.show()
    
    # Plot for the "decompress_one_sample" command
    plt.figure(figsize=(8, 5))
    for chromosome, group in decompress_one_sample_df.groupby('Chromosome'):
        plt.plot(group['Sample_Size'], group['Memory_Used_MB'], marker='o', color=color_map[chromosome], label=f'Chr {chromosome}')
    plt.plot(ALL_decompress_one_sample_df['Sample_Size'], ALL_decompress_one_sample_df['Memory_Used_MB'],marker='o', color='black', label="Mean")
    plt.xlabel('Sample Size')
    plt.ylabel('RAM (MB)')
    plt.xscale("log")
    titre = algo +' - RAM vs Sample Size for Decompress_one_sample Command'
    plt.title(titre)
    plt.legend(title='Chromosomes')#,bbox_to_anchor=(1.2, 1.1))
    plt.grid(True)
    plt.savefig("E:\VM_data\Git\VCF_compression_benchmarking\Images\\" + titre, bbox_inches='tight')
    plt.show()
    
    # Plot for the "random_access" command
    plt.figure(figsize=(8, 5))
    for chromosome, group in random_access_df.groupby('Chromosome'):
        plt.plot(group['Sample_Size'], group['RAM_per_MB'], marker='o', color=color_map[chromosome], label=f'Chr {chromosome}')
    plt.plot(ALL_random_access_df['Sample_Size'], ALL_random_access_df['RAM_per_MB'],marker='o', color='black', label="Mean")
    plt.xlabel('Sample Size')
    plt.ylabel('RAM per MB')
    plt.xscale("log")
    titre = algo +' - RAM per MB vs Sample Size for Random Access Command (bis)'
    plt.title(titre)
    plt.legend(title='Chromosomes')#,bbox_to_anchor=(1.2, 1.1))
    plt.grid(True)
    plt.savefig("E:\VM_data\Git\VCF_compression_benchmarking\Images\\" + titre, bbox_inches='tight')
    plt.show()
    
    # Plot for the "random_access_one_sample" command
    plt.figure(figsize=(8, 5))
    for chromosome, group in random_access_one_sample_df.groupby('Chromosome'):
        plt.plot(group['Sample_Size'], group['Memory_Used_MB'], marker='o', color=color_map[chromosome], label=f'Chr {chromosome}')
    plt.plot(ALL_random_access_one_sample_df['Sample_Size'], ALL_random_access_one_sample_df['Memory_Used_MB'],marker='o', color='black', label="Mean")
    plt.xlabel('Sample Size')
    plt.ylabel('RAM')
    plt.xscale("log")
    titre = algo +' - RAM vs Sample Size for Random access_one_sample Command'
    plt.title(titre)
    plt.legend(title='Chromosomes')#,bbox_to_anchor=(1.2, 1.1))
    plt.grid(True)
    plt.savefig("E:\VM_data\Git\VCF_compression_benchmarking\Images\\" + titre, bbox_inches='tight')
    plt.show()
