# Bash script for comparative analysis of bacterial genomes
#
# ⚠️ DO NOT execute this script entirely at once!
# Copy and paste individual command lines into the Linux terminal as needed.
# This file uses the .sh extension only to enable Bash syntax highlighting in text editors.
#
# Author: Marcus Vinicius Canário Viana
# Repository: https://github.com/canarioviana/workshop_amsterdam_umc
#

############################################################
## SUMMARY
############################################################

############################################################
## 0) Working directory
############################################################

# Connect to the server
ssh username@ipaddress

# Go to workshop directory
cd /mnt/4tb_1/workshop_umc

# Create assembly directory
mkdir -p username/comparative_analysis

# Go to assembly directory
cd username/comparative_analysis

############################################################
## 1) Input genomes 
############################################################

############################################################
## Input genomes from NCBI GenBank

# Leptospira interrogans from GenBank
# https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=173
# Filters:
#   Annotated by RefSeq
#   Assembly level complete
#   Exclude atypical genomes
#   Exclude MAGs
#   Exclude Genomes from large multi-isolate projects
#   Total of 62 genomes
# Download table

# Second filter:
# A serovar name in the column "Organism Name" and maximum of two strains per serovar.
# Total of 24 genomes

# Input file list
# Create the file 1_assembly_ids.tsv contaning accession numbers and strains names separated by tab, one per line and no column names

# Send the file to the directory "comparative_analysis"
# Remove Windows line brakes from the file
sed -i 's/\r//' 1_assembly_ids.tsv

############################################################
## Assembly files directory

# Create the genomes directory
mkdir 1_genomes

# Activate Conda environment
conda activate datasets
# Download compressed dehydrated directory
datasets download genome accession \
--assembly-version latest \
--include genome \
--dehydrated \
--inputfile 1_assembly_ids.tsv \
--filename genome_data.zip
# Unzip genome_data.zip
unzip genome_data.zip
# Rehydrate directory
datasets rehydrate --directory .
# Move assembly files to 1_genomes
mv ncbi_dataset/data/*/*.fna 1_genomes
# Remove sufix from file names and change file extension
cd 1_genomes
rename 's/(.*?_.*?)_.*/$1.fasta/' *.fna
cd ..
# Delete temporary files and directory
rm -r genome_data.zip ncbi_dataset README.md md5sum.txt
# Rename files according to the file 1_assembly_ids.tsv
parallel --colsep "\t" -a 1_assembly_ids.tsv mv 1_genomes/{1}.fasta 1_genomes/{2}.fasta

############################################################
## Private input genomes

# Put the private genome assembly files in fasta format the directory 1_genomes
# The files should have the extension .fsa and no special charactes in the name


############################################################
## 2) Taxonomy and assembly quality control
############################################################

############################################################
## GTDB-Tk (~1 minute)
# Requires 64GB of RAM if the species is not identified by the ANI screening step
# Activate Conda environment
conda activate gtdbtk
# Run the program
gtdbtk classify_wf \
--cpus $(nproc --ignore=1) \
--extension .fasta \
--mash_db /db/gtdbtk \
--genome_dir 1_genomes \
--out_dir 2_gtdbtk
# Deactivate Conda environment
conda deactivate
# Copy and rename output files
if [ -f "2_gtdbtk/gtdbtk.bac120.summary.tsv" ]; then
    cp 2_gtdbtk/gtdbtk.bac120.summary.tsv 2_gtdbtk_bacteria.tsv
fi
if [ -f "2_gtdbtk/gtdbtk.ar53.summary.tsv" ]; then
    cp 2_gtdbtk/gtdbtk.ar53.summary.tsv 2_gtdbtk_archaea.tsv
fi
# Compress the output directory
zip -r 2_gtdbtk.zip 2_gtdbtk
# Delete the output directory
rm -r 2_gtdbtk


############################################################
## TYGS
# https://tygs.dsmz.de/user_requests/new
# Send the files from 1_genomes


############################################################
## CheckM2 (~2 minutes)
# Activate Conda environment
conda activate checkm2
# Run the program
checkm2 predict \
--threads $(nproc --ignore=1) \
-x fasta \
--input 1_genomes \
--output-directory 2_checkm
# Deactivate Conda environment
conda deactivate
# Copy and rename the output file
cp 2_checkm/quality_report.tsv 2_checkm2.tsv
# Delete the output directory
# rm -r 2_checkm


############################################################
## GUNC (~7 minutes)
# Create an output directory
mkdir 2_gunc 2_gunc_temp
# Activate Conda environment
conda activate gunc
# Run the program
gunc run \
--threads $(nproc --ignore=1) \
--contig_taxonomy_output \
--file_suffix .fasta \
--input_dir 1_genomes \
--temp_dir 2_gunc_temp \
--out_dir 2_gunc
# Copy and rename the output file
cp 2_gunc/*maxCSS_level.tsv 2_gunc.tsv
# Plotting the data
# Create an output directory
mkdir 2_gunc/plot
# Loop through a list of files
for file in 2_gunc/diamond_output/*.out; do\
   gunc plot \
   --contig_display_num 0\
   --diamond_file $file \
   --out_dir 2_gunc/plot;\
done
# Deactivate Conda environment
conda deactivate
# Compress the output directory
zip -r 2_gunc.zip 2_gunc
# Delete the output directory
# rm -r 2_gunc 2_gunc_temp 


############################################################
## QUAST

# Activate Conda environment
conda activate quast
# Run the program
quast.py -m 0 -o 2_quast 1_genomes/*.fasta
# Deactivate Conda environment
conda deactivate
cp 2_quast/transposed_report.tsv 2_quast.tsv
# Compress the output directory
zip -r 2_quast.zip 2_quast
# Delete the output directory
rm -r 2_quast

############################################################
## Create list of genomes for downwstream analysis

# Copy the file 1_assembly_ids.tsv to 3_selected_genomes.tsv
cp 1_assembly_ids.tsv 3_selected_genomes.tsv
# Open the file in a text editor (nano)
nano 3_selected_genomes.tsv
# Delete the line of any genome that failed the quality control filters
# Save the file and exit
# Press CTRL + O to save the file. Confirm with Enter.
# Press CTRL + X to exit nano.


############################################################
## 3) Directory of selected genomes
############################################################

############################################################
## Create a directory for the selected genomes

# Create the directory
mkdir 3_selected_genomes
# Copy the selected files
parallel --colsep "\t" -a 3_selected_genomes.tsv cp 1_genomes/{2}.fasta 3_selected_genomes/
# Compress the directory 3_selected_genomes
zip -r 3_selected_genomes.zip 3_selected_genomes
# Compress the directory 1_genomes
zip -r 1_genomes.zip 1_genomes
# Delete the directory 1_genomes
rm -r 1_genomes

############################################################
## 4) Genome annotation (Prokka OR PGAP)
############################################################

# ############################################################
# ## Prokka
# # Local annotation.

# # Create directory for genome annotation
# mkdir 4_genome_annotation

# # Activate Conda environment
# conda activate prokka
# # Loop through a list of files
# for file in 3_selected_genomes/*.fasta; do
#     #Extract file name
#     filename=${file##*/}
#     #Extract sample name
#     sample=${filename%%.*}
#     prokka \
#     --cpus $(nproc --ignore=1) \
#     --addgenes \
#     --centre "" \
#     --outdir 4_genome_annotation/${sample} \
#     --prefix ${sample} \
#     $file
# done
# # Deactivate Conda environment
# conda deactivate

# # Generate Prokka annotation summary files
# echo -e "Sample\tGene\tCDS\trRNA\ttRNA\ttmRNA" > 4_genome_annotation.tsv
# find 11_genome_annotation -type f -name "*.txt" | while read -r file; do
#     # Extract file name
#     filename=${file##*/}
#     # Extract sample name
#     sample=${filename%%.*}
#     # Convert lines to tab separated table lines
#     cds=$(grep -i "^CDS:" "$file" | cut -d':' -f2 | tr -d ' ')
#     gene=$(grep -i "^gene:" "$file" | cut -d':' -f2 | tr -d ' ')
#     rrna=$(grep -i "^rRNA:" "$file" | cut -d':' -f2 | tr -d ' ')
#     trna=$(grep -i "^tRNA:" "$file" | cut -d':' -f2 | tr -d ' ')
#     tmrna=$(grep -i "^tmRNA:" "$file" | cut -d':' -f2 | tr -d ' ')
#     # Add line to the output file
#     echo -e "${sample}\t${gene}\t${cds}\t${rrna}\t${trna}\t${tmrna}" >> 11_genome_annotation.tsv
# done

# # If you already had a genome annotated using Prokka, put the annotation directory inside 4_genome_annotation
# # The annotation directory must have the same prefix as the annotation files. For example: 4_genome_annotation/samplename/samplename.fsa  

# # Compress the output directory
# zip -r 4_genome_annotation.zip 4_genome_annotation


############################################################
## PGAP (In case you do not want to use Prokka)
# Anotation already available from NCBI Refseq database
# NCBI selects the best genome assemblies from GenBank for the RefSeq database, so it has fewer genomes

# Create directory for genome annotation
mkdir 4_genome_annotation

# Activate Conda environment
conda activate datasets
# Download compressed and dehydrated directory
datasets download genome accession \
--assembly-version latest \
--include genome,protein,cds,gbff,gff3 \
--dehydrated \
--inputfile 3_selected_genomes.tsv \
--filename genome_data.zip
# Unzip genome_data.zip
unzip genome_data.zip
# Rehydrate directory
datasets rehydrate --directory .
# Go to the genomes directory
cd ncbi_dataset/data
# Rename files
for sample in *; do
    if [ -d "$sample" ]; then cd $sample
    mv protein.faa "$sample".faa
    mv cds_from_genomic.fna "$sample".ffn
    mv genomic.gbff "$sample".gbk
    mv genomic.gff "$sample".gff
    rename 's/(.*?_.*?)_.*/$1.fsa/' *.fna
    cd ..
    fi
done
# Deactivate Conda environment
conda deactivate

# Move genomes directories (GC*/) to ../../4_genome_annotation
mv GC*/ ../../4_genome_annotation
# Go to main directory
cd ../../
# Rename files according to the file 2_selected_genomes.tsv
parallel --colsep "\t" -a 2_selected_genomes.tsv mv 4_genome_annotation/{1} 4_genome_annotation/{2}
parallel --colsep "\t" -a 2_selected_genomes.tsv rename 's/{1}/{2}/' 4_genome_annotation/{2}/*
# Delete temporary files and directory
rm -r genome_data.zip ncbi_dataset README.md md5sum.txt

# If you already had a genome annotated using PGAP, put the annotation directory inside 4_genome_annotation
# The annotation directory must have the same prefix as the annotation files. For example: 4_genome_annotation/samplename/samplename.fsa  

# Compress the output directory
zip -r 4_genome_annotation.zip 4_genome_annotation


############################################################
## 5) Mobilome analysis
############################################################

############################################################
## MOB-suite (~5 minutes)

# Create an output directory
mkdir 5_mobsuite

# Activate Conda environment
conda activate mob_suite
# Loop through a list of files
for file in 4_genome_annotation/*/*.fsa; do
    # Extract sample name (PGAP annotation)
    sample=$(basename $file .fsa)
    # Run the program
    mob_recon -n $(nproc --ignore=1) \
    --infile $file \
    --outdir 5_mobsuite/${sample}_mobsuite
done
# Deactivate Conda environment
conda deactivate

# Merge contig_report.families.tsv files
if find 5_mobsuite/ -maxdepth 2 -type f -name "mobtyper_results.txt" | grep -q .; then
    for file in 5_mobsuite/*/contig_report.txt; do \
        cat $file;\
        echo;\
    done > 5_mobsuite/5_mobsuite_contig_report.families.tsv
fi
# Merge mobtyper_results.tsv files
if find 5_mobsuite/ -maxdepth 2 -type f -name "mobtyper_results.txt" | grep -q .; then
    for file in 5_mobsuite/*/mobtyper_results.txt; do \
        cat $file; \
        echo; \
    done > 5_mobsuite/5_mobsuite_mobtyper_results_all.tsv
fi
# Merge contig_report.families.tsv files
if find 5_mobsuite/ -maxdepth 2 -type f -name "mge.report.txt" | grep -q .; then
    for file in 5_mobsuite/*/mge.report.txt; do \
         cat $file; \
         echo; \
    done > 5_mobsuite/5_mobsuite_mge.report_all.tsv
fi
# Copy merged mobtyper result file to main directory
if ls 5_mobsuite/5_mobsuite_mobtyper_results_all.tsv > /dev/null 2>&1; then
    cp 5_mobsuite/5_mobsuite_mobtyper_results_all.tsv .
fi

# Rename MOB-suite fasta files
# Loop through a list of directories
for dir in 5_mobsuite/*/; do
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

# Compress the output directory
zip -r 5_mobsuite.zip 5_mobsuite

# A plasmid classified on the primary cluster AF922 is the chromosome 2.

############################################################
## 6) Phylogenetic analysis
############################################################

# ############################################################
# ## Phylogeny outgroup - Prokka annotation

# # Create a the tab separated text file "6_assembly_ids_outgroup.tsv", containing the assembly ID and species name, one genome per line.
# # Send the file to the working directory
# # In case of using Corynebacterium diphtheriae NCTC11397 as outgroup:
# # echo -e GCA_902809765.1"\t"NCTC11397 > 6_assembly_ids_outgroup.tsv

# # Create output directory
# mkdir -p 6_phylogeny_outgroup

# # Download genome assembly
# # Activate Conda environment
# conda activate datasets
# # Download compressed dehydrated directory
# datasets download genome accession \
# --assembly-version latest \
# --include genome \
# --dehydrated \
# --inputfile 6_assembly_ids_outgroup.tsv \
# --filename genome_data.zip
# # Unzip genome_data.zip
# unzip genome_data.zip
# # Rehydrate directory
# datasets rehydrate --directory .
# # Move assembly files to 1_genomes
# mv ncbi_dataset/data/*/*.fna 6_phylogeny_outgroup
# # Remove sufix from file names and change file extension
# cd 6_phylogeny_outgroup
# rename 's/(.*?_.*?)_.*/$1.fasta/' *.fna
# cd ..
# # Delete temporary files and directory
# rm -r genome_data.zip ncbi_dataset README.md md5sum.txt
# # Rename files according to the file 1_assembly_ids.tsv
# parallel --colsep "\t" -a 6_assembly_ids_outgroup.tsv mv 6_phylogeny_outgroup/{1}.fasta 6_phylogeny_outgroup/{2}.fasta

# # Annotate outgroup genome using Prokka
# # Activate Conda environment
# conda activate prokka
# # Loop through a list of files
# for file in 6_phylogeny_outgroup/*.fasta; do
#     #Extract file name
#     filename=${file##*/}
#     #Extract sample name
#     sample=${filename%%.*}
#     prokka \
#     --cpus $(nproc --ignore=1) \
#     --addgenes \
#     --centre "" \
#     --outdir 6_phylogeny_outgroup/${sample} \
#     --prefix ${sample} \
#     $file
# done
# # Deactivate Conda environment
# conda deactivate

# # Compress the output directory
# zip -r 6_phylogeny_outgroup.zip 6_phylogeny_outgroup

############################################################
## Phylogeny outgroup - PGAP annotation from RefSeq database

# Create a the tab separated text file "6_assembly_ids_outgroup.tsv", containing the assembly ID and species name, one genome per line.
# Send the file to the working directory

# Download outgroup Leptonema illini DSM 21528 (GCF_000243335.1)
# Create genome list
echo -e GCF_000243335.1"\t"DSM21528 > 6_assembly_ids_outgroup.tsv

# Create output directory
mkdir -p 6_phylogeny_outgroup

# Activate Conda environment
conda activate datasets
# Download compressed and dehydrated directory
datasets download genome accession \
--assembly-version latest \
--include genome,protein,cds,gbff,gff3 \
--dehydrated \
--inputfile 6_assembly_ids_outgroup.tsv \
--filename genome_data.zip
# Unzip genome_data.zip
unzip genome_data.zip
# Rehydrate directory
datasets rehydrate --directory .
# Deactivate Conda environment
conda deactivate
# Go to the genomes directory
cd ncbi_dataset/data
# Rename files
for sample in *; do
    if [ -d "$sample" ]; then cd $sample
    mv protein.faa "$sample".faa
    mv cds_from_genomic.fna "$sample".ffn
    mv genomic.gbff "$sample".gbk
    mv genomic.gff "$sample".gff
    rename 's/(.*?_.*?)_.*/$1.fsa/' *.fna
    cd ..
    fi
done
# Move genomes directories (GC*/) to ../../6_phylogeny_outgroup
mv GC*/ ../../6_phylogeny_outgroup
# Go to main directory
cd ../../
# Delete temporary files and directory
rm -r genome_data.zip ncbi_dataset README.md md5sum.txt
# Rename files according to the file 6_assembly_ids_outgroup.tsv
parallel --colsep "\t" -a 6_assembly_ids_outgroup.tsv mv 6_phylogeny_outgroup/{1} 6_phylogeny_outgroup/{2}
parallel --colsep "\t" -a 6_assembly_ids_outgroup.tsv rename 's/{1}/{2}/' 6_phylogeny_outgroup/{2}/*

# Compress the output directory
zip -r 6_phylogeny_outgroup.zip 6_phylogeny_outgroup

# ############################################################
# ## Phylogeny using PhyloPhlAn

# # Create output directory
# mkdir -p 6_phylogeny_phylophlan/faa
# # Go to faa directory
# cd 6_phylogeny_phylophlan/faa
# # Create symbolic links of faa files to save disk space
# ln -s ../../4_genome_annotation/*/*.faa .
# ln -s ../../6_phylogeny_outgroup/*/*.faa .
# # Go to PhyloPhlAn directory
# cd ..

# # Activate Conda environment
# conda activate phylophlan
# # Generate configuration files
# phylophlan_write_default_configs.sh
# # Run PhyloPhlAn with phylophlan database and configuration file supermatrix_aa.cfg
# phylophlan --nproc $(nproc --ignore=1) -d phylophlan --diversity low -f supermatrix_aa.cfg -i faa
# # The output file faa.tre is the tree generated by FastTree and has node support values
# # Run PhyloPhlAn with amphora2 database and configuration file supermatrix_aa.cfg 
# phylophlan --nproc $(nproc --ignore=1) -d amphora2 --diversity low -f supermatrix_aa.cfg -i faa
# # The output file faa.tre is the tree generated by FastTree and has node support values
# # Deactivate Conda environment
# conda deactivate

# # Go to main directory
# cd ..


############################################################
## Phylogeny using OrthoFinder (for small datasets)

# Create output directory
mkdir 6_phylogeny_orthofinder

# Create symbolic links of faa files to save disk space
ln -s $PWD/6_phylogeny_outgroup/*/*.faa 6_phylogeny_orthofinder/
ln -s $PWD/4_genome_annotation/*/*.faa 6_phylogeny_orthofinder/
# Go to OrthoFinder directory
cd 6_phylogeny_orthofinder

# Calculate and use the number open file descriptors, to avoid the error "Too many open files"
required_limit=$(ls -1 *.faa | wc -l | awk '{print $1 * $1}')
current_limit=$(ulimit -Sn)
if [ "$required_limit" -gt "$current_limit" ]; then
    echo "Increasing the number of open file descriptors to ${required_limit} to avoid the error 'Too many open files'"
    ulimit -n "$required_limit"
fi

# Activate Conda environment
conda activate orthofinder
#Run OrthoFinder using FastTree for tree inference (FASTER, ~17 minutes)
orthofinder -t $(nproc --ignore=1) -M msa -T fasttree -f . 2>&1 | tee orthofinder_fastree_stdout.log
# #Run OrthoFinder using IQTREE for tree inference (SLOWER)
# orthofinder -t $(nproc --ignore=1) -M msa -T iqtree -f . 2>&1 | tee orthofinder_iqtree_stdout.log
# Deactivate Conda environment
conda deactivate

#Go to maind directory
cd ..

# Compress the output directory
zip -r 6_phylogeny_orthofinder.zip \
6_phylogeny_orthofinder/OrthoFinder/Results_*/Species_Tree \
6_phylogeny_orthofinder/OrthoFinder/Results_*/MultipleSequenceAlignments/SpeciesTreeAlignment.fa \
6_phylogeny_orthofinder/orthofinder_*.log \
6_phylogeny_orthofinder/OrthoFinder/Results_*/Log.txt

# The number of single-copy genes used for the phylogeny is informed in the log files in 6_phylogeny_orthofinder/orthofinder_*.log 

# Visualize the tree online using iTOL
# Input file: 6_phylogeny_orthofinder/OrthoFinder/Results_*/Species_Tree/SpeciesTree_rooted.txt
https://itol.embl.de/upload.cgi
# Reroot the tree using DSM21258
# Branch lenghts: Ignore
# Advanced -> Branch metadata display -> Bootstrap metadata: Display
# Advanced -> Branch metadata display -> Display range: 0.7 to 1
# Datasets -> Create a dataset -> Type: Text label, Label: Serovar -> Create dataset
    # Paste the columns from the file tree_annotation.xlsx -> Update and close
    # Label size factor: 0.8x
# The strains of serovar Canicola (LJ178 and 782) are in different clades!

############################################################
# 7) Comparative analysis
############################################################

############################################################
## Pyani (Calculate ANI)

# Create an output directory
mkdir 7_pyani

# Go to directory
cd 7_pyani

# Create soft links for the genomes
ln -s ../4_genome_annotation/*/*.fsa .
rename s/fsa$/fasta/ *.fsa

# Activate Conda environment
conda activate pyani
# Run pyani
average_nucleotide_identity.py --gmethod seaborn -m ANIb -i . -o pyani -g
# Deactivate Conda environment
conda deactivate

# Go back to working directory
cd ..
# Compress the output directory
zip -r 7_pyani.zip 7_pyani

############################################################
## Genome-to-Denome Distace Calculator (GGDC) (calculate dDDH)
#http://ggdc.dsmz.de/

############################################################
## PPanGGOLiN (Pangenome)

# Generate list of .gbk input files from PGAP
find $PWD/4_genome_annotation/ -name '*.gbk' | parallel echo -e "{/.}'\t'{}" > 7_gbk_list_pangenome.tsv

# Activate environment
conda activate ppanggolin
# PPanGGOLiN workflow
ppanggolin workflow --cpu $(nproc --ignore=1) --anno 7_gbk_list_pangenome.tsv -o 7_ppanggolin -f
# Rarefaction (pangenome development)
ppanggolin rarefaction --cpu $(nproc --ignore=1) -p 7_ppanggolin/pangenome.h5 -o 7_ppanggolin/rarefaction
# Tile plot with dendrogram
ppanggolin draw -p 7_ppanggolin/pangenome.h5 --tile_plot --add_dendrogram -o 7_ppanggolin/tile_plot_dendrogram
# Tile plot with no cloud partition
ppanggolin draw -p 7_ppanggolin/pangenome.h5 --tile_plot --nocloud -o 7_ppanggolin/tile_plot_no_cloud 
# Info
ppanggolin info -p 7_ppanggolin/pangenome.h5 --content > 7_ppanggolin/content.txt
ppanggolin info -p 7_ppanggolin/pangenome.h5 --metadata > 7_ppanggolin/metadata.txt
ppanggolin info -p 7_ppanggolin/pangenome.h5 --parameters > 7_ppanggolin/parameters.txt
ppanggolin info -p 7_ppanggolin/pangenome.h5 --status > 7_ppanggolin/status.txt
# Stats
ppanggolin write_pangenome -p 7_ppanggolin/pangenome.h5 --stats -o 7_ppanggolin/stats
# Fasta - Nucleotide
ppanggolin fasta -p 7_ppanggolin/pangenome.h5 -o 7_ppanggolin/ -f --gene_families all
# Fasta - Protein
ppanggolin fasta -p 7_ppanggolin/pangenome.h5 -o 7_ppanggolin/ -f --prot_families all

# Pangenome annotation - Directory
pandir=7_ppanggolin
anndir="$pandir/annotation"
mkdir -p "$anndir"

# Pangenome annotation - EggNOG-Mapper
mkdir -p "$anndir"/eggnogmapper
conda activate eggnog-mapper
emapper.py \
    --cpu $(nproc --ignore=1) \
    -i "$pandir"/all_protein_families.faa \
    --output all_protein_families \
    --output_dir "$anndir"/eggnogmapper
conda deactivate
grep -v '^##' "$anndir"/eggnogmapper/all_protein_families.emapper.annotations \
    > "$anndir"/eggnogmapper/all_protein_families.emapper.annotations.tsv

# Pangenome annotation - COG classifier
conda activate cogclassifier
COGclassifier \
    -t $(nproc --ignore=1) \
    -i "$pandir"/all_protein_families.faa \
    -o "$anndir"/cogclassifier
conda deactivate

# Pangenome annotation - dbCAN
mkdir -p "$anndir"/dbcan
conda activate dbcan
run_dbcan CAZyme_annotation \
    --threads $(nproc --ignore=1) \
    --db_dir /db/dbcan \
    --mode protein \
    --input_raw_data "$pandir"/all_protein_families.faa \
    --output_dir "$anndir"/dbcan
conda deactivate

# Pangenome annotation - AMR Finder
mkdir -p "$anndir"/amrfinder
conda activate amrfinder
amrfinder \
    --threads $(nproc --ignore=1) \
    --database /db/amrfinder/latest \
    --plus \
    --annotation_format prokka \
    --protein "$pandir"/all_protein_families.faa \
    --protein_output "$anndir"/amrfinder/amrfinder.faa \
    --output "$anndir"/amrfinder/amrfinder.tsv
conda deactivate

# Pangenome annotation - RGI
mkdir -p "$anndir"/rgi
conda activate rgi
rgi main \
    -n $(nproc --ignore=1) \
    --clean \
    -t protein \
    -i "$pandir"/all_protein_families.faa \
    -o "$anndir"/rgi/rgi
conda deactivate

# Pangenome annotation - Diamond-VFDB
mkdir -p "$anndir"/vfdb
conda activate diamond
diamond \
    blastp \
    -p $(nproc --ignore=1) \
    --id 70 \
    --query-cover 70 \
    --evalue 1e-10 \
    --max-target-seqs 1 \
    -d /mnt/4tb_1/db/vfdb_full/vfdb_full.dmnd \
    -q "$pandir"/all_protein_families.faa \
    -o "$anndir"/vfdb/vfdb.tsv \
    -f 6 qseqid sseqid pident qcovhsp evalue bitscore mismatch gapopen qstart qend sstart send stitle
conda deactivate
awk -F'\t' \
'BEGIN{print "qseqid\tsseqid\tpident\tqcovhsp\tevalue\tbitscore\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tstitle"} \
{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13}' \
"$anndir"/vfdb/vfdb.tsv \
> "$anndir"/vfdb/vfdb_header.tsv

# Compress PPanGGLiN directory
zip -r "$pandir".zip "$pandir"