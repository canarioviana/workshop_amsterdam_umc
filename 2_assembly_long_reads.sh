#!/bin/bash
# Bash script for bacterial genome assembly from Oxford NanoPore sequencing data
#
# Author: Marcus Vinicius Canário Viana
# Date: 17/02/2025
# Repository: https://github.com/canarioviana/workshop_amsterdam_umc


############################################################
## SUMMARY OF END-TO-END GENOME ASSEMBLY WORKFLOW FROM LONG-READS
############################################################

## 0) Working directory
## 1) Sequencing reads directory and files
    # Reads stored as local files
## 2) Raw reads quality assessment
    # NanoPlot
## 3) Raw reads trimming
    # Fastplong
## 4) Trimmed reads quality assessment
    # NanoPlot
## 5) De novo assembly
    # Flye
## 6) Organization of de novo assembly files
## 7) Assembly quality assessment
    # CheckM2
    # GUNC
    # QUAST
    # Calculation of vertical sequencing coverage
## 8) Taxonomic assignment
    # GTDB-Tk
    # TYGS
    # JspeciesWS
## 9) Plasmid identification
    # MOB-suite
## 10) Assignment of contigs to molecules
    # MOB-suite and an inhouse script


############################################################
## 0) Working directory
############################################################

# Connect to the server
ssh user@ipaddress

# Go to workshop directory
cd /mnt/4tb_1/workshop_umc

# Create assembly directory
mkdir -p username/long_reads

# Go to assembly directory
cd username/long_reads


############################################################
## 1) Sequencing reads directory and files
############################################################

############################################################
## Reads from NCBI SRA database

# Sequencing data to download
# Leptospira interrogans strain N116
# BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA828002
# BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN27628310
# Genome: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023515975.1/
# SRA: https://www.ncbi.nlm.nih.gov/sra/SRX27382703
# NanoPore flow cell: R10.3
# Basecalling model: Super accuracy model (SUP, performed with Guppy (v5.0.11+2b6dbff))
# Raw reads accession: SRR32032837

# Declare accession id
accession=SRR32032837
# Declare sample name
sample=N116

# Create output directory
mkdir 1_reads

# SRA Tools
# Activate Conda environment
conda activate sra-tools

# Run prefetch
prefetch -p -O 1_reads "${accession}"

# Run fasterq-dump
cd 1_reads
fasterq-dump \
--threads $(nproc --ignore=1) \
-p \
--split-files "${accession}" \
-O .

# Deactivate Conda environment
conda deactivate

# Delete temporary directory
rm -r "${accession}"

# Compress files
pigz -p $(nproc --ignore=1) "${accession}"*.fastq

# Rename files
mv "${accession}.fastq.gz" "${sample}.fq.gz"

# Back to assembly directory
cd ..


############################################################
## 2) Raw reads quality assessment
############################################################

############################################################
## NanoPlot

# Create an output directory
mkdir 2_nanoplot

# Activate Conda environment
conda activate nanoplot

# Loop through a list of files
for reads in 1_reads/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz

    # Verify if the input files are empty
    if [ ! -s "$reads" ]; then
        echo "Warning: The input files of sample ${sample} are empty. Skiping..."
        echo -e "${sample}" >> 2_nanoplot_skiped_samples.tsv
        # Skip to the next iteration
        continue
    fi

    echo "NanoPlot is processing sample: ${sample}"
    mkdir "2_nanoplot/${sample}"
    NanoPlot -t $(nproc --ignore=1) --fastq ${reads} -p "${sample}_" -o "2_nanoplot/${sample}"
done

# Deactivate Conda environment
conda deactivate

# Compress the output directory
zip -r 2_nanoplot.zip 2_nanoplot

# Delete the output directory
rm -r 2_nanoplot


############################################################
## 3) Raw reads trimming
############################################################

############################################################
## Fastplong

# Create an output directory
mkdir 3_fastplong

# Activate Conda environment
conda activate fastplong
# Loop through a list of files
for reads in 1_reads/*.fq.gz; do

    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz

    # Inform the current sample being processed
    echo "Processing sample: ${sample}"
    echo "Reads file ${reads}"

    # Run Fastplong
    fastplong \
    --thread $(nproc --ignore=1) \
    --disable_quality_filtering \
    --length_required 500 \
    -i "$reads" \
    -o 3_fastplong/"${sample}_trimmed.fq.gz" \
    --html 3_fastplong/"$sample"_fastp.html \
    --json 3_fastplong/"$sample"_fastp.json
done

# Deactivate Conda environment
conda deactivate

# Compress report files
# zip -r 3_fastplong.zip 3_fastplong/*.json 3_fastplong/*.html
# Delete report files
# rm 3_fastplong/*.json 3_fastplong/*.html


############################################################
## 4) Trimmed reads quality assessment
############################################################

############################################################
## NanoPlot

# Create an output directory
mkdir 4_nanoplot

# Activate Conda environment
conda activate nanoplot

# Declare input directory
input_directory="3_fastplong"
# Loop through a list of files
for reads in "${input_directory}"/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz

    # Verify if the input files are empty
    if [ ! -s "$reads" ]; then
        echo "Warning: The input files of sample ${sample} are empty. Skiping..."
        echo -e "${sample}" >> 4_nanoplot_skiped_samples.tsv
        # Skip to the next iteration
        continue
    fi

    echo "NanoPlot is processing sample: ${sample}"
    echo "A warning is expected if the Conda environment was installed by the root."
    mkdir "4_nanoplot/${sample}"
    NanoPlot -t $(nproc --ignore=1) --fastq ${reads} -p "${sample}_" -o "4_nanoplot/${sample}"
done

# Deactivate Conda environment
conda deactivate

# Compress the output directory
zip -r 4_nanoplot.zip 4_nanoplot

# Delete the output directory
rm -r 4_nanoplot


############################################################
## 5) De novo assembly
############################################################

############################################################
## Flye (~37 minutes)

# Create an output directory
mkdir 5_flye

# Activate Conda environment
conda activate flye

# Set library type
# valid_options: --pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw | --nano-corr | --nano-hq
seq_type="--nano-hq"

# Declare input directory
input_directory="3_fastplong"
# Loop through a list of files
for reads in "${input_directory}"/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz
    # Run Flye
    echo "Assembling the genome of sample: $sample"
    flye -t $(nproc --ignore=1) \
    "${seq_type}" \
    "$reads" \
    -o 5_flye/"$sample"_flye
done

# Deactivate Conda environment
conda deactivate

# Compress output files
zip -r 5_flye.zip 5_flye

############################################################
## Flye -> Medaka

# Create an output directory
mkdir 5_flye_medaka

# Activate Conda environment
conda activate medaka

# Avoid using the GPU
CUDA_VISIBLE_DEVICES="" 

# Declare inference model according to your sequencing library
# List the available inference models using the command: medaka tools list_models
# Current sample has NanoPore flow cell R10.3 and basecalling model SUP
model=r103_sup_g507

# Loop through a list of files
for reads in 3_fastplong/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz

    # Run medaka
    medaka_consensus \
    -i ${reads} \
    -d 5_flye/"${sample}"_flye/assembly.fasta \
    -o 5_flye_medaka/"${sample}"_medaka \
    -t $(nproc --ignore=1) \
    -m "${model}"
done


############################################################
## 6) Organization of de novo assembly files
############################################################

############################################################
## Create assemblies directory

mkdir 6_assemblies

############################################################
## Organizing assemblies from Flye -> Medaka

# Loop through a list of directories
for dir in 5_flye_medaka/*/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_medaka*}
    # Copy and rename the assembly file
    cp "${dir}consensus.fasta" "6_assemblies/${sample}.fasta"
done

############################################################
## Compress the assemblies directory

zip -r 6_assemblies.zip 6_assemblies

############################################################
## Delete the assemblers output directories

rm -r 5_flye


############################################################
## 7) Assembly quality assessment
############################################################

############################################################
## CheckM2 (~1 minute)

# Activate Conda environment
conda activate checkm2

# Run the program
checkm2 predict \
--threads $(nproc --ignore=1) \
-x fasta \
--input 6_assemblies \
--output-directory 7_checkm

# Deactivate Conda environment
conda deactivate

# Copy and rename the output file
cp 7_checkm/quality_report.tsv 7_checkm2.tsv

# Delete the output directory
rm -r 7_checkm

############################################################
## GUNC (~2 minutes)

# Create an output directory
mkdir 7_gunc 7_gunc_temp

# Activate Conda environment
conda activate gunc

# Run the program
gunc run \
--threads $(nproc --ignore=1) \
--contig_taxonomy_output \
--file_suffix .fasta \
--input_dir 6_assemblies \
--temp_dir 7_gunc_temp \
--out_dir 7_gunc

# Copy and rename the output file
cp 7_gunc/*maxCSS_level.tsv 7_gunc.tsv

# Plotting the data
# Create an output directory
mkdir 7_gunc/plot

# Loop through a list of files
for file in 7_gunc/diamond_output/*.out; do\
   gunc plot \
   --contig_display_num 0\
   --diamond_file $file \
   --out_dir 7_gunc/plot;\
done

# Deactivate Conda environment
conda deactivate

# Compress the output directory
zip -r 7_gunc.zip 7_gunc

# Delete the output directory
rm -r 7_gunc 7_gunc_temp 

############################################################
## QUAST

# Activate Conda environment
conda activate quast

# Run the program
quast.py -m 0 -o 7_quast 6_assemblies/*.fasta

# Deactivate Conda environment
conda deactivate
cp 7_quast/transposed_report.tsv 7_quast.tsv

# Compress the output directory
zip -r 7_quast.zip 7_quast

# Delete the output directory
rm -r 7_quast

############################################################
## Vertical sequencing coverage

# Declare input directory
reads_directory="3_fastplong"
# Create output file
echo -e Sample"\t"Coverage > 7_coverage.tsv
# Loop through a list of files
for file in 6_assemblies/*.fasta; do
    # Extract assembly name
    assembly=$(basename $file .fasta)
    # Extract assembly name
    sample=${assembly%%_*}
    # Inform the sample
    echo -e Calculating sequencing coverage for assembly: $assembly
    echo -e Assembly file: ${file}
    echo -e Read file: "${reads_directory}"/${sample}_trimmed.fq.gz
    # Count bases in assembly
    bases_in_assembly=$(grep -v '^>' "$file" | tr -d '\n' | wc -c)
    echo -e Bases in assembly: $bases_in_assembly
    # Count bases in sequencing files
    # bases_in_reads=$(zcat "${reads_directory}"/"$sample"*.gz | awk 'NR%4==2 {print $0}' | tr -d '\n' | wc -c)
    bases_in_reads=$(zcat "${reads_directory}"/${sample}_trimmed.fq.gz | 
                     awk 'NR%4==2 {print length}' | paste -sd+ | bc)
    echo -e Bases in reads: $bases_in_reads
    # Calculate vertical coverage
    if [ "$bases_in_assembly" -gt 0 ]; then
        # coverage=$(($bases_in_reads / $bases_in_assembly))
        coverage=$(echo "scale=2; $bases_in_reads / $bases_in_assembly" | bc)
        echo -e Coverage: "${coverage}\n"
        echo -e "${assembly}\t${coverage}" >> 7_coverage.tsv
    else
        echo -e "${assembly}\t0" >> 7_coverage.tsv
    fi
done


############################################################
## 8) Taxonomic assignment
############################################################

############################################################
## GTDB-Tk

# Requires 64GB of RAM if the species is not identified by the ANI screening step
# Activate Conda environment
conda activate gtdbtk

# Run the program
gtdbtk classify_wf \
--cpus $(nproc --ignore=1) \
--extension .fasta \
--mash_db /db/gtdbtk \
--genome_dir 6_assemblies \
--out_dir 8_gtdbtk

# Deactivate Conda environment
conda deactivate

# Copy and rename output files
if [ -f "8_gtdbtk/gtdbtk.bac120.summary.tsv" ]; then
    cp 8_gtdbtk/gtdbtk.bac120.summary.tsv 8_gtdbtk_bacteria.tsv
fi
if [ -f "8_gtdbtk/gtdbtk.ar53.summary.tsv" ]; then
    cp 8_gtdbtk/gtdbtk.ar53.summary.tsv 8_gtdbtk_archaea.tsv
fi

# Compress the output directory
zip -r 8_gtdbtk.zip 8_gtdbtk

# Delete the output directory
rm -r 8_gtdbtk


############################################################
## TYGS

# https://tygs.dsmz.de/user_requests/new
# Send the files from 6_assemblies


############################################################
## JspeciesWS

# Tetra Correlation Search (TCS)
# https://jspecies.ribohost.com/jspeciesws/#analyse

# Session code for previous result: F6E4DCCF5CDC31D53C1C

############################################################
## 9) Plasmid identification
############################################################

############################################################
## MOB-suite

# Create an output directory
mkdir 9_mobsuite

# Activate Conda environment
conda activate mob_suite

# Loop through a list of files
for file in 6_assemblies/*.fasta; do
    # Extract sample name
    sample=$(basename $file .fasta)
    # Run the program
    mob_recon -n $(nproc --ignore=1) \
    --infile $file \
    --outdir 9_mobsuite/${sample}_mobsuite
done

# Deactivate Conda environment
conda deactivate

# Merge contig_report.txt files
# Delete merged file if it exists
if [ -f 9_mobsuite/contig_report_all.tsv ]; then
    rm 9_mobsuite/contig_report_all.tsv
fi
# Initialize control variable to check in the header was printed
header_printed=0
# Check if any contig_report.txt files exist
if find 9_mobsuite/ -maxdepth 2 -type f -name "contig_report.txt" | grep -q .; then
    # Iterate over all found contig_report.txt files
    for file in 9_mobsuite/*/contig_report.txt; do
        # Test if the header was not printed yet
        if [ $header_printed -eq 0 ]; then
            # If not, contatenate the entire file
            cat "$file" >> 9_mobsuite/contig_report_all.tsv
            # Mark that the header was printed
            header_printed=1
        else
        # The header was printed, so contanate file and ignore its header
        tail -n +2 "$file" >> 9_mobsuite/contig_report_all.tsv  
        fi
        # Add a new line to separate the results of each sample
        # echo >> 9_mobsuite/contig_report_all.tsv
    done
fi

# Merge mobtyper_results.txt files
# Delete merged file if it exists
if [ -f 9_mobsuite/mobtyper_results.txt ]; then
    rm 9_mobsuite/mobtyper_results_all.tsv
fi
# Initialize control variable to check in the header was printed
header_printed=0
# Check if any mobtyper_results.txt files exist
if find 9_mobsuite/ -maxdepth 2 -type f -name "mobtyper_results.txt" | grep -q .; then
    # Iterate over all found mobtyper_results.txt files
    for file in 9_mobsuite/*/mobtyper_results.txt; do
        # Test if the header was not printed yet
        if [ $header_printed -eq 0 ]; then
            # If not, contatenate the entire file
            cat "$file" >> 9_mobsuite/mobtyper_results_all.tsv
            # Mark that the header was printed
            header_printed=1
        else
        # The header was printed, so contanate file and ignore its header
        tail -n +2 "$file" >> 9_mobsuite/mobtyper_results_all.tsv  
        fi
        # Add a new line to separate the results of each sample
        # echo >> 9_mobsuite/mobtyper_results_all.tsv
    done
fi

# Merge mge.report.txt files
# Delete merged file if it exists
if [ -f 9_mobsuite/mge.report_all.tsv ]; then
    rm 9_mobsuite/mge.report_all.tsv
fi
# Initialize control variable to check in the header was printed
header_printed=0
# Check if any mge.report.txt files exist
if find 9_mobsuite/ -maxdepth 2 -type f -name "mge.report.txt" | grep -q .; then
    # Iterate over all found mge.report.txt files
    for file in 9_mobsuite/*/mge.report.txt; do
        # Test if the header was not printed yet
        if [ $header_printed -eq 0 ]; then
            # If not, contatenate the entire file
            cat "$file" >> 9_mobsuite/mge.report_all.tsv
            # Mark that the header was printed
            header_printed=1
        else
        # The header was printed, so contanate file and ignore its header
        tail -n +2 "$file" >> 9_mobsuite/mge.report_all.tsv  
        fi
        # Add a new line to separate the results of each sample
        # echo >> 9_mobsuite/mge.report_all.tsv
    done
fi

# Copy merged mobtyper result file to main directory
if ls 9_mobsuite/mobtyper_results_all.tsv > /dev/null 2>&1; then
    cp 9_mobsuite/mobtyper_results_all.tsv 9_mobsuite_mobtyper_results_all.tsv
fi

# Rename MOB-suite fasta files
# Loop through a list of directories
for dir in 9_mobsuite/*/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_mobsuite*}
    # Go to input directory
    # Add sample name to the chromosome file name 
    mv ${dir}chromosome.fasta ${dir}"$sample"_chromosome.fasta
    # Add sample name to the plasmid file name 
    rename s/plasmid_/"$sample"_plasmid_/ ${dir}plasmid*
    # Leave input diretory
done


#########################################################################
## 10) Assignment of contigs to molecules
############################################################

############################################################
## Add molecule attribution to contigs. Required for batch genome submission.

# Copy assemblies directory
cp -r 6_assemblies 10_assemblies_for_analysis
# Change extension to .fsa
rename "s/.fasta$/.fsa/" 10_assemblies_for_analysis/*.fasta
# Loop through a list of files (assemblies)
for assembly in 10_assemblies_for_analysis/*.fsa; do
    # Extract file name
    filename=${assembly##*/}
    # Extract sample name
    sample=${filename%%.*}
    # Inform sample
    echo "Analyzing $sample"
    # Loop through a list of files (mob-suite result files)
    awk '1' 9_mobsuite/"$sample"_mobsuite/contig_report.txt | while IFS=$'\t' read -r sample_id molecule_type primary_cluster_id secondary_cluster_id contig_id others; do
        # Ignore column names
        if [ "$sample_id" == "sample_id" ]; then
            continue
        fi
        # Conditional: plasmid contig
        if [ "$molecule_type" == "plasmid" ]; then
            # Create new contig header
            new_contig_id=$(echo "$contig_id" "[plasmid-name="p"$primary_cluster_id""]")
            # Inform sample and contig
            echo Sample: "$sample" - Contig: "$new_contig_id"
            # Replace old header with the new one
            sed -i "s/^>${contig_id}/>${new_contig_id}/" $assembly
        fi
        # Conditional: chromosome contig
        if [ "$molecule_type" == "chromosome" ]; then
            # Create new contig header
            new_contig_id=$(echo "$contig_id" "[chromosome=1]")
            # Inform sample and contig
            echo Sample: "$sample" - Contig: "$new_contig_id"
            # Replace old header with the new one
            sed -i "s/^>${contig_id}/>${new_contig_id}/" $assembly
        fi
    done
done
# Compress the output directory
zip -r 10_assemblies_for_analysis.zip 10_assemblies_for_analysis
# Compress mob-suite directory
zip -r 9_mobsuite.zip 9_mobsuite
# Delete output file
rm -r 9_mobsuite

# Rename molecules from strain N116
sed -i 's/plasmid-name=pAF922/chromosome=2/' 10_assemblies_for_analysis/N116.fsa
sed -i 's/pAA693/pl_URK_4a/' 10_assemblies_for_analysis/N116.fsa
sed -i 's/pAA694/pl_URK_4b/' 10_assemblies_for_analysis/N116.fsa
sed -i 's/pAA696/pl_URK_4c/' 10_assemblies_for_analysis/N116.fsa

# In case of a novel plasmid, you will have to change its temporary name given by MOB-suite to an appropriate and shorter name.


############################################################
# Back to the main user directory

cd /mnt/4tb_1/workshop_umc/username


############################################################
## Genome submission to GenBank

# Submission portal
# https://submit.ncbi.nlm.nih.gov/subs/

# Create BioProject
# https://submit.ncbi.nlm.nih.gov/subs/bioproject/

# Create Biosample
# https://submit.ncbi.nlm.nih.gov/subs/biosample/
# BioSample template (Required form)
# https://submit.ncbi.nlm.nih.gov/biosample/template/

# Submit the genome
# https://submit.ncbi.nlm.nih.gov/subs/genome/
# Choose the "Batch submission", even for a single submission: “New submission” escolher a opção “Batch/multiple genomes (maximum 400 per submission)”.
# Choose the genome annotation with PGAP: "Annotate this prokaryotic genome in the NCBI Prokaryotic Annotation Pipeline (PGAP) before its release"
# Upload a table containing the genomes metadata. A template (Batch genomes: Genome Info file template) is available in https://submit.ncbi.nlm.nih.gov/templates/.
# Upload the .fsa files in the directory 10_assemblies_for_analysis

# Submit sequencing reads
# https://submit.ncbi.nlm.nih.gov/subs/sra/
# Upload a table containing the sequencing metadata. A template (SRA: Metadata spreadsheet with sample names) is available in https://submit.ncbi.nlm.nih.gov/templates/.
# Upload the .fq.gz files in directory 1_reads

