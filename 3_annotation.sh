# Bash script for bacterial genome annotation
#
# ⚠️ DO NOT execute this script entirely at once!
# Copy and paste individual command lines into the Linux terminal as needed.
# This file uses the .sh extension only to enable Bash syntax highlighting in text editors.
#
# Author: Marcus Vinicius Canário Viana
# Date: 17/02/2025
# Repository: https://github.com/canarioviana/workshop_amsterdam_umc/
# More info: see README.md in the repository


############################################################
## GENOME ANNOTATION WORKFLOW
############################################################
## This script is the continuation of the script "bacterial_genome_assembly.sh"


############################################################
## 11) Genome annotation
############################################################

#############################################################
## Genome annotation with Prokka

# The assembly files should be in the directory 10_assemblies_for_analysis with the sufix .fsa and no special characters in their names.

# Activate Conda environment
conda activate prokka
# Loop through a list of files
for file in 10_assemblies_for_analysis/*.fsa; do
    #Extract file name
    filename=${file##*/}
    #Extract sample name
    prefix=${filename%%.*}
    prokka \
    --cpus $(nproc --ignore=1) \
    --addgenes \
    --centre "" \
    --outdir 11_genome_annotation/${prefix} \
    --prefix ${prefix} \
    $file
done
# Deactivate Conda environment
conda deactivate

# Generate annotation summary files
echo -e "Sample\tGene\tCDS\trRNA\ttRNA\ttmRNA" > 11_genome_annotation.tsv
find 11_genome_annotation -type f -name "*.txt" | while read -r file; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%.*}
    # Convert lines to tab separated table lines
    cds=$(grep -i "^CDS:" "$file" | cut -d':' -f2 | tr -d ' ')
    gene=$(grep -i "^gene:" "$file" | cut -d':' -f2 | tr -d ' ')
    rrna=$(grep -i "^rRNA:" "$file" | cut -d':' -f2 | tr -d ' ')
    trna=$(grep -i "^tRNA:" "$file" | cut -d':' -f2 | tr -d ' ')
    tmrna=$(grep -i "^tmRNA:" "$file" | cut -d':' -f2 | tr -d ' ')
    # Add line to the output file
    echo -e "${sample}\t${gene}\t${cds}\t${rrna}\t${trna}\t${tmrna}" >> 11_genome_annotation.tsv
done

# Compress the output directory
zip -r 11_genome_annotation.zip 11_genome_annotation
# Delete the output directory
# rm -r 11_genome_annotation


############################################################
## 12) Genome annotation - Virulence and resistance genes
############################################################

############################################################
## AMRFinderPlus

# Create an output directory
mkdir -p 12_amrfinder

# Activate Conda environment
conda activate amrfinder
# Loop through a list of sample files
for file in 11_genome_annotation/*/*.faa; do
    # Extract sample name
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%.faa}
    
    # Inform current sample
    echo "AMRFinderPlus is processing sample: ${sample}"

    # Method 1 - Run AMRFinderPlus using nucleotide and proteins sequences
    # Run main software
    amrfinder \
    --threads $(nproc --ignore=1) \
    --database /db/amrfinder/latest \
    --plus \
    --annotation_format prokka \
    --nucleotide "11_genome_annotation/${sample}/${sample}.fsa" \
    --protein "11_genome_annotation/${sample}/${sample}.faa"  \
    --gff "11_genome_annotation/${sample}/${sample}.gff" \
    --name "${sample}" \
    --nucleotide_output "12_amrfinder/${sample}_amrfinder.fasta"\
    --protein_output "12_amrfinder/${sample}_amrfinder.faa"\
    --output "12_amrfinder/${sample}_amrfinder.tsv"

    # # Method 2 - Run AMRFinderPlus using only proteins sequences
    # # Run main software
    # amrfinder \
    # --threads $(nproc --ignore=1) \
    # --database /db/amrfinder/latest \
    # --plus \
    # --protein "${file}"  \
    # --name "${sample}" \
    # --protein_output "12_amrfinder/${sample}_amrfinder.faa"\
    # --output "12_amrfinder/${sample}_amrfinder.tsv"
done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
zip -r 12_amrfinder.zip 12_amrfinder
# Delete the output directory
rm -r 12_amrfinder

############################################################
## ResFinder (Antimicrobial resistance prediction)
# http://genepi.food.dtu.dk/resfinder

############################################################
## RGI (Antimicrobial resistance prediction)
# https://card.mcmaster.ca/analyze/rgi

# Activate Conda environment
conda activate rgi
mkdir 12_rgi
# Create file to save the results table
echo -e "Genome\tORF_ID\tContig\tStart\tStop\tOrientation\tCut_Off\tPass_Bitscore\tBest_Hit_Bitscore\tBest_Hit_ARO\tBest_Identities\tARO\tModel_type\tSNPs_in_Best_Hit_ARO\tOther_SNPs\tDrug Class\tResistance Mechanism\tAMR Gene Family\tPredicted_DNA\tPredicted_Protein\tCARD_Protein_Sequence\tPercentage Length of Reference Sequence\tID\tModel_ID\tNudged\tNote\tHit_Start\tHit_End\tAntibiotic" > 12_rgi/rgi_all.txt
# Loop through a list of files
for file in 11_genome_annotation/*/*.faa; do
    sample=$(basename ${file} .faa)
    #Run RGI
    rgi main -n $(nproc --ignore=1) --clean -t protein -i $file -o 12_rgi/"$sample"
    sed -i "s/^/$sample\t/" 12_rgi/"$sample".txt #Sample name
    tail -n +2 12_rgi/"$sample".txt >> 12_rgi/rgi_all.txt #Results
    echo >> 12_rgi/rgi_all.txt #Empty line
done
# Deactivate Conda environment
conda deactivate

# Copy the summary file to main directory
cp 12_rgi/rgi_all.txt 12_rgi.tsv
# Compress the output directory
zip -r 12_rgi.zip 12_rgi
# Delete the output directory
rm -r 12_rgi

############################################################
## VFanalyzer (Virulence gene prediction)
# http://www.mgc.ac.cn/VFs/


############################################################
## 12) Genome annotation - Mobile Genetic Elements
############################################################

############################################################
## CRISPRcasFinder local (CRISPR/Cas system prediction)

# Create the output directory
mkdir 12_crisprcasfinder
# Go to the output directory
cd 12_crisprcasfinder
# Copy binaries to the current directory
cp -r /db/CRISPRCasFinder/* .

# Activate Conda environment
conda activate crisprcasfinder
# Install dependencies
macsydata install -u CASFinder==3.1.0 # You only have to run this command the first time you use this software
# Loop through a list of files
for file in ../11_genome_annotation/*/*.fsa; do
    prefix=$(basename ${file} .fasta)
    # Run CRISPR-CasFinder
    perl CRISPRCasFinder.pl -cas -html -in $file -out "${prefix}"
done
# Deactivate Conda environment
conda deactivate

# Remove the CRISPRCasFinder files from the current directory
rm -r CasFinder-2.0.3 ccf.environment.yml COPYING COPYRIGHT CRISPRCasFinder* crispr-cas_logo.png *.sh *.txt install_test README.md sel392v2.so singularity supplementary_files supplementary_tools

# Concatenate CRISPR-CasFinder summaries
> CRISPR-Cas_summary_all.tsv
for directory in */; do
    sample=$(basename $directory /)
    echo -e $sample >> CRISPR-Cas_summary_all.tsv
    cat "$sample"/TSV/CRISPR-Cas_summary.tsv >> CRISPR-Cas_summary_all.tsv
    echo -e "\n" >> CRISPR-Cas_summary_all.tsv
done

# Edit the summary file header
sed -i 's/Evidence-levels/Nb_arrays_evidence-level_1\tNb_arrays_evidence-level_2\tNb_arrays_evidence-level_3\tNb_arrays_evidence-level_4/g' CRISPR-Cas_summary_all.tsv
# Separate into columns
sed -i 's/,Nb/\tNb/g' CRISPR-Cas_summary_all.tsv
# Remove the prefix
sed -i 's/Nb_arrays_evidence-level_.=//g' CRISPR-Cas_summary_all.tsv
# Remove commas
sed -i 's/,\t/\t/g' CRISPR-Cas_summary_all.tsv
# Leave the directory
cd ..

# Copy the summary file to main directory
cp 12_crisprcasfinder/CRISPR-Cas_summary_all.tsv 12_crisprcasfinder_summary.tsv
# Compress the output directory
zip -r 12_crisprcasfinder.zip 12_crisprcasfinder
# Delete the output directory
rm -r 12_crisprcasfinder

############################################################
## CRISPRcasFinder online (CRISPR/Cas system prediction)
# https://crisprcas.i2bc.paris-saclay.fr/CrisprCasFinder/Index
# Send the chromosome file (*.fasta)

############################################################
## PHASTER (Prophage prediction)

# https://phaster.ca/
# Send the GenBank file (*.gb, *.gbk or *.gbff) or chromosome file (*.fasta)

############################################################
## PHASTEST (Prophage prediction)

# https://phastest.ca/
# Send the GenBank file (*.gb, *.gbk or *.gbff) or the chromosome file (*.fasta)

############################################################
## IslandViewer
# https://www.pathogenomics.sfu.ca/islandviewer/upload/


############################################################
## 12) Genome annotation - Molecular typing
############################################################

############################################################
## FastMLST (Multi-Locus Sequence Typing)
# https://github.com/EnzoAndree/FastMLST

# Update the database using the command below. It requires the same user that created the Conda environment.
# fastmlst -t 1 --update-mlst

# Create  ourput directory
mkdir 12_fastmlst

# Activate Conda environment
conda activate fastmlst
# Run the software
fastmlst \
-t 1 \
-s '\t' \
-to 12_fastmlst/fastmlst.tsv \
-fo 12_fastmlst/mlst_concat_alleles.fasta \
--splited-output 12_fastmlst/mlst_alleles \
11_genome_annotation/*/*.fsa
# Deactivate Conda environment
conda deactivate

# Remove the "." from the sample name
sed -i s/.fsa// 12_fastmlst/fastmlst.tsv
# Copy the summary file to main directory
cp 12_fastmlst/fastmlst.tsv 12_fastmlst.tsv
# Compress the output directory
zip -r 12_fastmlst.zip 12_fastmlst
# Delete the output directory
rm -r 12_fastmlst

############################################################
## PubMLST

# https://pubmlst.org/organisms
# Enter the organism page, enter Typing page
# "Query a sequence" -> "Single sequence"
# Send the chromosome file (*.fasta) and click on "Submit"


############################################################
## 12) Genome annotation - Functional annotation
############################################################

############################################################
## EggNOG-mapper (Gene functional annotation)

# Create the output directory
mkdir -p 12_eggnog_mapper

# Activate Conda environment
conda activate eggnog-mapper
# Loop through a list of files
for file in 11_genome_annotation/*/*.faa; do
    prefix=$(basename ${file} .faa)
    # Run EggNOG-mapper
    emapper.py --cpu $(nproc --ignore=1) -i $file --output "${prefix}" --output_dir 12_eggnog_mapper
    tail +5 12_eggnog_mapper/"${prefix}".emapper.annotations > 12_eggnog_mapper/"${prefix}".emapper.annotations.tsv
done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
zip -r 12_eggnog_mapper.zip 12_eggnog_mapper
# Delete the output directory
rm -r 12_eggnog_mapper emappertmp*

############################################################
## COG classifier (Gene functional annotation)

# Create an output directory
mkdir 12_cogclassifier

# Activate Conda environment
conda activate cogclassifier
# Loop through a list of files
for file in 11_genome_annotation/*/*.faa; do
    prefix=$(basename ${file} .faa)
    # Run COG classifier
    COGclassifier -t $(nproc --ignore=1) -i $file -o 12_cogclassifier/"${prefix}"
done
# Deactivate Conda environment
conda deactivate

# Compress the output directory
zip -r 12_cogclassifier.zip 12_cogclassifier
# Delete the output directory
rm -r 12_cogclassifier

############################################################
## Blastkoala (Gene functional annotation)

# https://www.kegg.jp/blastkoala/

############################################################
## dbCAN 

# Create an output directory
mkdir -p 12_dbcan

# Activate Conda environment
conda activate dbcan
# Loop through a list of sample files
for file in 11_genome_annotation/*/*.faa; do
    # Extract sample name
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%.faa}

    # Create output directory
    mkdir -p "12_dbcan/${sample}_dbcan"

    # Run main software
    run_dbcan CAZyme_annotation \
    --threads $(nproc --ignore=1) \
    --db_dir /db/dbcan \
    --mode protein \
    --input_raw_data $file \
    --output_dir "12_dbcan/${sample}_dbcan"

    # Rename output files
    for outputfile in 12_dbcan/${sample}_dbcan/*; do
            outputfilename=${outputfile##*/}
        mv ${outputfile} "12_dbcan/${sample}_dbcan/${sample}_${outputfilename}"
    done
done
# Deactivate Conda environment
conda deactivate

# Compress directory
echo "Compressing output directory"
zip -r 112_dbcan.zip 12_dbcan

############################################################
## 13) Genome annotation - Circular map
############################################################

# https://proksee.ca/