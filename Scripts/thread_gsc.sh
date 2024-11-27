#!/bin/bash
echo "Start"

# Chromosome range
CHROMOSOMES=(1 22 10 15 5 1 22 10 15 5) # 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2)

# Sample sizes to crop
sample_sizes=(2500 2250 2000 1500 1000 900 800 750 700 600 500 250 200 150 100 50 25 15 10 5 2 1) #

#Number of thread to use
thread_list=(1 2 3 4 8)

# CSV file for logging results
output_time_RAM_csv="gsc_usage_thread_time_RAM.csv"
output_filesize_csv="gsc_usage_thread_filesize.csv"
echo "Command,Chromosome,Sample_Size,Threads,Time_Taken,Memory_Used_MB,Input_file_size" > "$output_time_RAM_csv"
echo "Command,Chromosome,Sample_Size,Threads,Input_file_size,Output_file_size" > "$output_filesize_csv"


# Iterate over each chromosome
for CHR in "${CHROMOSOMES[@]}"; do

    DIRECTORY="../1000gp_p3/ALL.chr${CHR}"
    mkdir -p "$DIRECTORY"
    CROPPED_VCF="$DIRECTORY/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf"
    gzip -d -k -c "../1000gp_p3/zip_files/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" > "$CROPPED_VCF"
    echo "gzip done"
    
    for size in "${sample_sizes[@]}"; do
        INTERMEDIAIRE_VCF="$DIRECTORY/ALL.chr${CHR}.${size}_samples.vcf"
        
        sample_list=$(../tools/bcftools/bcftools query -l "$CROPPED_VCF" | shuf | head -n "$size" | tr '\n' ',' | sed 's/,$//')
	 echo "CROPPED VCF == "
	 echo "$CROPPED_VCF"

        # Crop and filter VCF files
        bcftools view -s "$sample_list" -o "$INTERMEDIAIRE_VCF" -O v "$CROPPED_VCF"
	 rm "$CROPPED_VCF"

	 CROPPED_VCF="$DIRECTORY/chr${CHR}.${size}_samples.vcf"
        bcftools view -i'GT="alt"' -Ov -o "$CROPPED_VCF" "$INTERMEDIAIRE_VCF"
		
        rm "$INTERMEDIAIRE_VCF"
	 for thread in "${thread_list[@]}"; do

		for i in {1..3}; do
		    SIZE_CROPPED_VCF=$(stat -c%s "$CROPPED_VCF")

		    # Run GSC Compress
		    COMPRESSED_VCF="$DIRECTORY/chr${CHR}.${size}_samples.vcf.gsc"

		    # Run the command and capture the output in a variable
		    /usr/bin/time --output="$output_time_RAM_csv" -a -f "compress,$CHR,$size,$thread,%E,%M,$SIZE_CROPPED_VCF" ./../tools/GSC/gsc compress --in "$CROPPED_VCF" --out "$COMPRESSED_VCF" --threads "$thread"
		    
		    SIZE_COMPRESSED_VCF=$(stat -c%s "$COMPRESSED_VCF")
		    echo "compress,$CHR,$size,$thread,$SIZE_CROPPED_VCF,$SIZE_COMPRESSED_VCF" >> "$output_filesize_csv"

		    rm "$COMPRESSED_VCF"
		    rm -f "$COMPRESSED_VCF.tmp"
		done
		echo "Processed the $size samples on chromosome $CHR with $thread thread(s)"
	 done
    done
    
	 rm "$CROPPED_VCF"

done


