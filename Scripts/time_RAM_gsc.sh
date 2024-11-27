#!/bin/bash
echo "Start"

# Chromosome range
CHROMOSOMES=(1 22 1 22 10 10 15 15) # 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2)

# Sample sizes to crop
sample_sizes=(2500 2250 2000 1500 1000 900 800 750 700 600 500 250 200 150 100 50 25 15 10 5 2 1)

# CSV file for logging results
output_csv="gsc_usage_log_chrom_G.csv"
echo "Command,Chromosome,Sample_Size,Time_Taken,Memory_Used_MB,Input_file_size" > "$output_csv"

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
	
	for i in {1..3}; do
		SIZE_CROPPED_VCF=$(stat -c%s "$CROPPED_VCF")

		# Run GSC Compress
		COMPRESSED_VCF="$DIRECTORY/chr${CHR}.${size}_samples.vcf.gsc"

		# Run the command and capture the output in a variable
		/usr/bin/time --output="$output_csv" -a -f "compress,$CHR,$size,%E,%M,$SIZE_CROPPED_VCF" ./../tools/GSC/gsc compress --in "$CROPPED_VCF" --out "$COMPRESSED_VCF"
		
		SIZE_COMPRESSED_VCF=$(stat -c%s "$COMPRESSED_VCF")
		
		# Run GSC Decompress
		DECOMPRESSED_VCF="$DIRECTORY/chr${CHR}.${size}_samples.decompressed.vcf"
		/usr/bin/time --output="$output_csv" -a -f "decompress,$CHR,$size,%E,%M,$SIZE_COMPRESSED_VCF" ./../tools/GSC/gsc decompress --in "$COMPRESSED_VCF" --out "$DECOMPRESSED_VCF"
		
		# Random Access for 2 random positions
		positions=($(awk '!/^#/ { 
		    if (count < 3) { 
		        pos[count] = $2 
		    } else {
		        idx = int(rand() * (count + 1))
		        if (idx < 3) {
		            pos[idx] = $2
		        }
		    }
		    count++
		} END { for (i = 0; i < 3; i++) print pos[i] }' "$CROPPED_VCF"))
		one_sample=$(../tools/bcftools/bcftools query -l "$CROPPED_VCF" | shuf | head -n 1 | tr '\n' ',' | sed 's/,$//')
		echo "SAMPLE"
		echo "$one_sample"
		
		for pos in "${positions[@]}"; do
		    echo "$pos"
		    RANDOM_ACCESS_VCF="$DIRECTORY/chr${CHR}.${size}_samples.decompressed.${pos}.vcf"
		    /usr/bin/time --output="$output_csv" -a -f "random_access,$CHR,$size,%E,%M,$SIZE_COMPRESSED_VCF" ./../tools/GSC/gsc decompress -M --range ${CHR}:${pos},${pos} --in "$COMPRESSED_VCF" --out "$RANDOM_ACCESS_VCF"
		    
		    ONE_SAMPLE_RANDOM_ACCESS_VCF="$DIRECTORY/chr${CHR}.${size}_samples.decompressed.${pos}_${one_sample}.vcf"
		    /usr/bin/time --output="$output_csv" -a -f "random_access_one_sample,$CHR,$size,%E,%M,$SIZE_COMPRESSED_VCF" ./../tools/GSC/gsc decompress -M --range ${CHR}:${pos},${pos} -s "$one_sample" --in "$COMPRESSED_VCF" --out "$ONE_SAMPLE_RANDOM_ACCESS_VCF"
		    
		    rm "$RANDOM_ACCESS_VCF"
		    rm "$ONE_SAMPLE_RANDOM_ACCESS_VCF"
		done

		ONE_SAMPLE_DECOMPRESSION_VCF="$DIRECTORY/chr${CHR}.${size}_samples.decompressed.${one_sample}.vcf"
		/usr/bin/time --output="$output_csv" -a -f "decompress_one_sample,$CHR,$size,%E,%M,$SIZE_COMPRESSED_VCF" ./../tools/GSC/gsc decompress -M -s "$one_sample" --in "$COMPRESSED_VCF" --out "$ONE_SAMPLE_DECOMPRESSION_VCF"
		    
		
		rm "$COMPRESSED_VCF"
		rm -f "$COMPRESSED_VCF.tmp"
		rm "$DECOMPRESSED_VCF"
		rm "$ONE_SAMPLE_DECOMPRESSION"
	done

        echo "Processed the $size samples on chromosome $CHR"
    done
    
	rm "$CROPPED_VCF"
done
