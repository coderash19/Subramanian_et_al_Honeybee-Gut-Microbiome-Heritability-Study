# ====================================================================================
# PATRILINE ASSESSMENT 
# ====================================================================================

# Step 1. Used Windows GUI to generate .Dat file as per 'COLONY' software manual and placed it in the COLONY Linux package folder 

# Step 2. Used Linux to run .Dat file from inside the 'COLONY' Linux package folder for honey bee colony 1
cd /home/Colony/colony2_Lnx_25_02_2024 # location of 'COLONY' software; run all commands from inside the 'COLONY' folder.
./colony2s.ifort.out IFN:COLONY-C1-rerun.Dat # input the .dat file name after IFN:

# Repeated above steps for my honey bee colony two samples with COLONY-C2.Dat file


# ====================================================================================
# MICROBIOME - PROCESSING PIPELINE
# ====================================================================================

# This study was conducted with two honey bee colonies. For readability, the code below references only Colony One in folder/file names.
# Please note that this pipeline was run for Colony One samples first, and then the same code was run for Colony Two samples. 

# ----------------------------------------------------------------------------------------------------------------------------------------------------
# PART 1. BASECALLING with Dorado
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# This folder needs to be in the directory 
dna_r10.4.1_e8.2_400bps_sup@v4.3.0


# First, create a bash script for the SLURM job with the following information: text header, module load, and basecaller loop code for all file names in the folder
# ----------------------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash
#SBATCH --account=def-zayed
#SBATCH --gpus-per-node=v100l:2
#SBATCH --mem=32G
#SBATCH --time=00-12:00         # time (DD-HH:MM)
#SBATCH --mail-user=as19@yorku.ca
#SBATCH --mail-type=ALL

module load dorado/0.6.1

filename='Filenames.txt' #contains list of all 
exec 4<$filename
echo Start
while read -u4 p ; do
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.3.0 /path/to/colony_one/pod5/$p --emit-fastq --device cuda:all > $p.colony_one_basecall.fastq
done


# Then, run the bash script from the directory where dna_r10.4.1_e8.2_400bps_sup@v4.3.0 folder is
# ----------------------------------------------------------------------------------------------------------------------------------------------------

cd /path/to/colony_one 
sbatch "${BASH SCRIPT NAME}".sh

# After job completion, transferred fastq files from DAC server to lab server to a folder called Colony_one_fastq


# Combine all fastq files together
# ----------------------------------------------------------------------------------------------------------------------------------------------------

cd /path/to/Fastq/Colony_one_fastq #cd to where I moved the basecall files in the server
cat *.fastq > Combined_Colony_One.fastq # combine them all to a new fastq file in same folder

# ----------------------------------------------------------------------------------------------------------------------------------------------------
# PART 2. DEMULTIPLEXING  
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Demultiplex with Minibar python - no trimming yet (ran in home directory in minibar folder) 
# File containing primer and index info formatted according to minibar expectations: C1-indexprimer-info/C1-minibar-demux-info.txt

cd /path/to/minibar-master # folder that contains Minibar loaded

python minibar.py C1-indexprimer-info/C1-minibar-demux-info.txt /path/to/Fastq/Colony_one_fastq/Combined_Colony_One.fastq -e 4 -E 4 -l 80 -F -S

# NOTE: 
# index edit distance (-e) 4
# primer edit distance (-E) 4


# ----------------------------------------------------------------------------------------------------------------------------------------------------
# PART 3. ADAPTER AND PRIMER TRIMMING
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Step 1: trim primers from 5' end using -g option  
# ------------------------------------------------------------------------------------------------------------

# Prepared .fasta file that lists PRIMER sequences to be trimmed
# Moved demultiplexed files to /path/to/Fastq/C1-demultiplexed folder
# Created .txt file that contain the sample id for each sample in the colony listed one per row

filename='/path/to/Fastq/C1-demultiplexed/sampleids-C1-rerun.txt'
exec 4<$filename
echo Start
while read -u4 p ; do
/usr/local/bin/cutadapt /path/to/Fastq/C1-demultiplexed/$p.fastq -j 10 -e 4 -g file:/path/to/important-files/PRIMERS.fasta -o /path/to/Fastq/C1-demultiplexed/C1-trimmed-5/$p.{name}.5.fastq
done

# move unknown files to unknown directory
source_dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-5"
dest_dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-5/C1-trimmed-5unk"
for file in $(find "$source_dir" -type f -name "*unknown*"); do
    mv "$file" "$dest_dir"
done

# cat F and R files for each 5' --> 3' trimmed sample 
sample_file="/path/to/Fastq/C1-demultiplexed/sampleids-C1-rerun.txt"
source_dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-5"
dest_dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-5/C1-trimmed-5FR"

while IFS= read -r sample; do
    file1="${source_dir}/${sample}.Forward-F515.5.fastq"
    file2="${source_dir}/${sample}.Reverse-R806.5.fastq"
    output_file="${dest_dir}/${sample}_5FR.fastq"
    cat "$file1" "$file2" > "$output_file"
done < "$sample_file"


# Step 2: trim primers from the 3' end using -a option
# ------------------------------------------------------------------------------------------------------------

# Same steps as above using the trimmed files from the 5' end

filename='/path/to/Fastq/C1-demultiplexed/C1-trimmed-5/C1-trimmed-5FR/sampleids-C1-rerun-5FR.txt'
exec 4<$filename
echo Start
while read -u4 p ; do
/usr/local/bin/cutadapt /data4/AishwaryaProject/Fastq/C1-demultiplexed/C1-trimmed-5/C1-trimmed-5FR/$p.fastq -j 10 -e 4 -a file:/data4/AishwaryaProject/important-files/RC-PRIMERS.fasta -o /data4/AishwaryaProject/Fastq/C1-demultiplexed/C1-trimmed-53/$p.{name}.53.fastq
done

# move unknown files to unknown directory
source_dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-53"
dest_dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-53/C1-trimmed-53unk"
for file in $(find "$source_dir" -type f -name "*unknown*"); do
    mv "$file" "$dest_dir"
done

# cat F and R files for each 5'-3' trimmed sample 
sample_file="/path/to/Fastq/C1-demultiplexed/C1-trimmed-5/C1-trimmed-5FR/sampleids-C1-rerun-5FR.txt"
source_dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-53"
dest_dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-53/C1-trimmed-53FR"

while IFS= read -r sample; do
    file1="${source_dir}/${sample}.RC-Forward-F515.53.fastq"
    file2="${source_dir}/${sample}.RC-Reverse-R806.53.fastq"
    output_file="${dest_dir}/${sample}_53FR.fastq"
    cat "$file1" "$file2" > "$output_file"
done < "$sample_file"


# Step 3: Rename files in the 53FR folder (intermediate step for a more convenient file naming convention)
# ------------------------------------------------------------------------------------------------------------

# remove _5FR from name
dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-53/C1-trimmed-53FR"
for file in "$dir"/*_5FR*; do
    new_file="${file//_5FR/}"
    mv "$file" "$new_file"
done

# replace _53FR with final
dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-53/C1-trimmed-53FR"
for file in "$dir"/*_53FR*; do
    new_file="${file//_53FR/_final}"
    mv "$file" "$new_file"
done


# ----------------------------------------------------------------------------------------------------------------------------------------------------
# PART 4. SIZE FILTERING with Cutadapt 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Size filter +/- 10 base pairs for all trimmed files

dir="/path/to/Fastq/C1-demultiplexed/C1-trimmed-53/C1-trimmed-53FR"

for file in "$dir"/*.fastq; do
    base_name="${file%.fastq}"
    output_file="${base_name}_filter.fastq"
    /usr/local/bin/cutadapt "$file" -m 240 -M 260 -o "$output_file"
done


# ----------------------------------------------------------------------------------------------------------------------------------------------------
# PART 5. INSTALLING AND PREPPING FOR with qiime2
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Install in home directory on cedar
# ----------------------------------------------------------------------------------------------------------------------------------------------------

cd /path/to/home

module load apptainer

apptainer build qiime2-2021.11.sif docker://quay.io/qiime2/core:2021.11


# Rename all the files to the  .+_.+_L[0-9][0-9][0-9]_R[12]_001\.fastq\.gz format that is required by qiime2
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# I did it one at a time but can be done with a bash script...

# .+_.+_L[0-9][0-9][0-9]_R[12]_001\.fastq\.gz
# G1_1_L000_R1_001.fastq.gz
# sample_G1_100_final_filter.fastq.gz

# ----------------------------------------------------------------------------------------------------------------------------------------------------
# PART 6. QIIME2 Steps
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Step 1. qiime2 tools import - convert fasta.gz to qza 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --time=00-24:00         # time (DD-HH:MM)
#SBATCH --mem=64G
#SBATCH --mail-user=as19@yorku.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-zayed

module load apptainer

# Path to the directory containing fastq.gz files
INPUT_DIR="/path/to/colony_one/final-filter-colony-one"

# Path to the temporary directory
TEMP_DIR="/path/to/colony_one/final-filter-colony-one/temporary_import"

# Create the temporary directory
mkdir -p "${TEMP_DIR}"

# Loop through each fastq.gz file in the input directory
for file in "${INPUT_DIR}"/*.fastq.gz; do
    # Extract the filename without extension
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.fastq.gz}"  # Remove both .fastq.gz

    # Create a subdirectory in the temporary directory for each sample
    sample_id="${filename_no_ext}"
    mkdir -p "${TEMP_DIR}/${sample_id}"
    
    # Copy the fastq.gz file into the subdirectory
    cp "${file}" "${TEMP_DIR}/${sample_id}/${filename}"
done

# Path to the directory where QIIME 2 artifacts will be saved
OUTPUT_DIR="/path/to/colony_one/qiime2-tools-import"

# Import each file into QIIME 2 and specify output based on filename
for subdir in "${TEMP_DIR}"/*; do
    if [ -d "$subdir" ]; then
        # Get the sample ID from the directory name
        sample_id=$(basename -- "$subdir")
        
        # Import the files into QIIME 2 using the temporary directory structure
        apptainer exec qiime2-2021.11.sif qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path "$subdir" \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt \
            --output-path "${OUTPUT_DIR}/${sample_id}.qza"
    fi
done

# Remove the temporary directory
rm -r "${TEMP_DIR}"


# Step 2. qiime2 vsearch dereplicate-sequences
# ----------------------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --time=00-24:00         # time (DD-HH:MM)
#SBATCH --mem=64G
#SBATCH --mail-user=as19@yorku.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-zayed

module load apptainer

# Directory where qza files are located
input_dir="/path/to/colony_one/qiime2-tools-import"
output_dir="/path/to/colony_one/qiime2-vsearch-dereplicated-seqs"

# Loop through all qza files in the input directory
for file in ${input_dir}/*.qza; do
    if [[ -f "$file" ]]; then
        filename=$(basename -- "$file")
        
        # Remove "_L000_R1_001" from filename
        filename_noext="${filename%_L000_R1_001.qza}.qza"

        # Construct output filenames
        output_table="${output_dir}/${filename_noext}-table.qza"
        output_rep_seqs="${output_dir}/${filename_noext}-rep-seqs.qza"

        # Run the QIIME command
        apptainer exec qiime2-2021.11.sif qiime vsearch dereplicate-sequences \
          --i-sequences "$file" \
          --o-dereplicated-table "$output_table" \
          --o-dereplicated-sequences "$output_rep_seqs"
    fi
done


# Step 3. qiime2 cluster-features-closed-reference 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --time=00-30:00         # time (DD-HH:MM)
#SBATCH --mem=64G
#SBATCH --mail-user=as19@yorku.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-zayed

module load apptainer

# Define the base directory path
base_dir="/path/to/colony_one/qiime2-vsearch-dereplicated-seqs"

# Loop through each rep-seqs.qza file in the base directory
for rep_seq_file in "$base_dir"/*-rep-seqs.qza; do
    # Extract the sample name by removing the suffix
    sample=$(basename "$rep_seq_file" -rep-seqs.qza)
    
    # Define the corresponding table file
    table_file="$base_dir/${sample}-table.qza"
    
    # Define output directory for clustered results
    output_dir="/path/to/colony_one/qiime2-cluster-features-closed-ref/qiime2-silva-cluster-features-closed-ref-c1"
    
    # Run the QIIME 2 command using the apptainer
    apptainer exec qiime2-2021.11.sif qiime vsearch cluster-features-closed-reference \
        --i-table "$table_file" \
        --i-sequences "$rep_seq_file" \
        --i-reference-sequences /scratch/as19/silva-138-99-seqs-515-806.qza \
        --p-perc-identity 0.98 \
        --p-strand 'both' \
        --o-clustered-table "$output_dir/${sample}-clustered-table-98-silva.qza" \
        --o-clustered-sequences "$output_dir/${sample}-clustered-sequences-98-silva.qza" \
        --o-unmatched-sequences "$output_dir/${sample}-unmatched-sequences-98-silva.qza"
done


# 4. qiime2 feature-classifier classify-consensus-vsearch 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --time=00-30:00         # time (DD-HH:MM)
#SBATCH --mem=64G
#SBATCH --mail-user=as19@yorku.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-zayed

module load apptainer

# Define the base directory path
base_dir="/path/to/colony_one/qiime2-cluster-features-closed-ref/qiime2-silva-cluster-features-closed-ref-c1"

# Define output directory for classification results
output_dir="/path/to/colony_one/qiime2-feature-classify-cons-vsearch/qiime2-silva-feature-classify-cons-vsearch-c1"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each clustered-sequences.qza file in the base directory
for clustered_sequences_file in "$base_dir"/*-clustered-sequences-98-silva.qza; do
    # Extract the sample name by removing the suffix
    sample=$(basename "$clustered_sequences_file" -clustered-sequences-98-silva.qza)

    # Run the QIIME feature-classifier classify-consensus-vsearch command
    apptainer exec qiime2-2021.11.sif qiime feature-classifier classify-consensus-vsearch \
      --i-query "$clustered_sequences_file" \
      --i-reference-reads /scratch/as19/silva-138-99-seqs-515-806.qza \
      --i-reference-taxonomy /scratch/as19/silva-138-99-tax-515-806.qza \
      --p-perc-identity 0.98 \
      --p-strand 'both' \
      --p-threads 10 \
      --o-classification "$output_dir/${sample}-classification-98-silva.qza"
done


# 5. qiime2 taxa collapse (family)
# ----------------------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --time=00-24:00         # time (DD-HH:MM)
#SBATCH --mem=64G
#SBATCH --mail-user=as19@yorku.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-zayed

module load apptainer

# Define the base directory path
base_dir="/path/to/colony_one/qiime2-feature-classify-cons-vsearch/qiime2-silva-feature-classify-cons-vsearch-c1"

# Define output directory for classification results
output_dir="/path/to/colony_one/qiime2-taxa-collapse/qiime2-silva-taxa-collapse-family-c1"

# Loop through each clustered-sequences.qza file in the base directory
for classify_file in "$base_dir"/*-classification-98-silva.qza; do
    # Extract the sample name by removing the suffix
    sample=$(basename "$classify_file" -classification-98-silva.qza)

    apptainer exec qiime2-2021.11.sif qiime taxa collapse \
        --i-table "/path/to/colony_one/qiime2-cluster-features-closed-ref/qiime2-silva-cluster-features-closed-ref-c1/${sample}-clustered-table-98-silva.qza" \
        --i-taxonomy "$classify_file" \
        --p-level 5 \
        --o-collapsed-table "$output_dir/${sample}-family-silva.qza"
done


# Step 6. qiime2 taxa collapse (genus)
# ----------------------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --time=00-24:00         # time (DD-HH:MM)
#SBATCH --mem=64G
#SBATCH --mail-user=as19@yorku.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-zayed

module load apptainer

# Define the base directory path
base_dir="/path/to/colony_one/qiime2-feature-classify-cons-vsearch/qiime2-silva-feature-classify-cons-vsearch-c1"

# Define output directory for classification results
output_dir="/path/to/colony_one/qiime2-taxa-collapse/qiime2-silva-taxa-collapse-genus-c1"

# Loop through each clustered-sequences.qza file in the base directory
for classify_file in "$base_dir"/*-classification-98-silva.qza; do
    # Extract the sample name by removing the suffix
    sample=$(basename "$classify_file" -classification-98-silva.qza)

    apptainer exec qiime2-2021.11.sif qiime taxa collapse \
        --i-table "/path/to/colony_one/qiime2-cluster-features-closed-ref/qiime2-silva-cluster-features-closed-ref-c1/${sample}-clustered-table-98-silva.qza" \
        --i-taxonomy "$classify_file" \
        --p-level 6 \
        --o-collapsed-table "$output_dir/${sample}-genus-silva.qza"
done

# ====================================================================================
# NEXT STEPS 
# ====================================================================================

# Next steps are conducted on R. This includes: (1) converting .qza files to .csv files for R, (2) Mistag filtering, and (3) Statistical Analysis and Figure Plotting

# Please see the file "Subramanian_et_al_R_ANALYSES.sh" for the following R pipeline used.

