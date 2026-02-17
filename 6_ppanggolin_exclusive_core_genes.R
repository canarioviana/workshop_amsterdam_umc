# # Install required packages
# install.packages("dplyr")
# install.packages("tidyverse")

# ===============================
# Libraries
# ===============================
library(dplyr)
library(tidyr)

# ===============================
# Output directory
# ===============================

# Create R Studio output directory
dir.create("exclusive_genes")

# ===============================
# Load Data
# ===============================

# Create the file genome_groups.tsv with the columns "Genome" for genome ID and "Group" for group name. 
# Put the files genome_groups.tsv and matrix_annotated.tsv in the working directory

# Import gene groups (matrix.csv)
matrix <- read.delim(
  file = "matrix.csv",
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
# Count the number of genome columns
first_genome_col <- 15 # Last column before the genomes columns.
num_genomes <- ncol(matrix) - first_genome_col + 1
rm(first_genome_col)
rm(matrix)

# Import gene groups (matrix_annotated.tsv)
matrix_annotated <- read.delim(
  file = "matrix_annotated.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Interval of genome columns from matrix_annotated.tsv
first_genome_col <- 15
last_genome_col <- first_genome_col + num_genomes - 1
total_cols <- ncol(matrix_annotated)

# Create a dataframe with gene presence table
gene_groups <- matrix_annotated %>% 
  select(1, all_of(first_genome_col):all_of(last_genome_col))

# Create a dataframe with gene annotation
gene_annotation <- matrix_annotated %>% 
  select(1, all_of(last_genome_col + 1):all_of(total_cols))

# Delete data to free memory
rm(matrix_annotated)
gc()

# Import genome groups (genome_groups.tsv)
genome_groups <- read.delim(
  file = "genome_groups.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Dataframe with the Genome and Group columns
genome_groups <- genome_groups %>% 
  select(Genome, Group) %>% 
  mutate(Group = if_else(is.na(Group) | Group == "", "NA_group", as.character(Group)))

# ===============================
# Long format
# ===============================
gene_long <- gene_groups %>%
  pivot_longer(
    cols = -1,
    names_to = "Genome",
    values_to = "Genes"
  ) %>%
  mutate(Present = !is.na(Genes) & Genes != "")
# "Genome", "Genes" and "Present" are column names for gene_long

# Delete data to free memory
rm(gene_groups)
gc()

# ===============================
# Add group info
# ===============================
gene_long <- gene_long %>%
  left_join(genome_groups, by = "Genome")

# ===============================
# Count genomes per group
# ===============================
group_sizes <- genome_groups %>%
  count(Group, name = "n_genomes")

# ===============================
# Presence stats per Gene & Group
# ===============================
gene_group_stats <- gene_long %>%
  group_by(Gene, Group) %>%
  summarise(
    n_present = sum(Present),
    .groups = "drop"
  ) %>%
  left_join(group_sizes, by = "Group") %>%
  mutate(pct_present = n_present / n_genomes)

# ===============================
# Frequency of core genes (Decimal, from 0 to 1)
# ===============================
min_core_pct <- 1

# ===============================
# Identify core-exclusive genes
# ===============================
core_exclusive_stats <- gene_group_stats %>%
  group_by(Gene) %>%
  filter(
    # Condition 1: minimum percentage in the target group
    any(pct_present >= min_core_pct) &
      # Condition 2: no presence in all other groups
      sum(n_present > 0) == 1
  ) %>%
  ungroup() %>%
  filter(pct_present >= min_core_pct)

# Add annotation
core_exclusive_stats_annotated <- core_exclusive_stats %>%
  left_join(gene_annotation, by = "Gene")

# Write table
write.table(
  core_exclusive_stats_annotated,
  file = "exclusive_genes/core_exclusive_stats_annotated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Target gene list
target_genes_list <- core_exclusive_stats %>% select(Gene, Group)

# Delete data to free the memory
rm(group_sizes, gene_group_stats, core_exclusive_stats)
gc()

# ===============================
# Extract genes
# ===============================
core_exclusive_genes <- gene_long %>%
  inner_join(target_genes_list, by = c("Gene", "Group")) %>%
  filter(Present) %>%
  separate_rows(Genes, sep = ",\\s*")

# Add annotation
core_exclusive_genes_annotated <- core_exclusive_genes %>%
  left_join(gene_annotation, by = "Gene")

# Delete data to free the memory
rm(gene_long, gene_annotation, target_genes_list)
gc()

# Write table
write.table(
  core_exclusive_genes_annotated,
  file = paste0(
    "exclusive_genes/core_exclusive_min",
    min_core_pct * 100,
    "pct_annotated.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Write table. One file per group
for (g in unique(core_exclusive_genes_annotated$Group)) {
  out <- core_exclusive_genes_annotated %>% filter(Group == g)
  
  write.table(
    out,
    file = paste0(
      "exclusive_genes/core_exclusive_",
      g,
      "_min",
      min_core_pct * 100,
      "pct_annotated.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
rm(out)

# ===============================
# Summary
# ===============================
core_summary_counts <- core_exclusive_genes %>% 
  distinct(Gene, Group) %>% 
  count(Group, name = "n_core_exclusive")

# Delete data to free memory
rm(core_exclusive_genes,genome_groups)
gc()

# Save summary
write.table(
  core_summary_counts,
  file = "exclusive_genes/exclusive_core_counts_per_group.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Delete data to free memory
rm(core_summary_counts, core_exclusive_genes_annotated)
gc()

# ===============================
# Get fasta sequences
# ===============================

# Extract exclusive genes ids
write.table(
  unique(core_exclusive_stats_annotated$Gene),
  file = "exclusive_genes/exclusive_genes_ids.tsv",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Delete data to free memory
rm(core_exclusive_stats_annotated)

# # Get exclusive genes sequences
# ssh username@ipaddress
# cd /mnt/4tb_1/workshop_umc/backup_plan/comparative_analysis/7_ppanggolin
# # Send the directory "exclusive_genes" to 7_ppanggolin
# # Extract the fasta files
# conda activate seqkit
# seqkit grep -f exclusive_genes/exclusive_genes_ids.tsv -r ../4_genome_annotation/*/*.faa > exclusive_genes/exclusive_genes.faa
# seqkit grep -f exclusive_genes/exclusive_genes_ids.tsv -r ../4_genome_annotation/*/*.ffn > exclusive_genes/exclusive_genes.ffn


