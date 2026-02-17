# # Install required packages
install.packages("tidyverse")

# Set the PPanGGOLiN output directory as Working Directory (Session -> Set Working Directory -> Choose Directory)
# From the PPanGGOLiN output you will only need the files "matrix.csv" and the directory "annotation" generated using the script "leptospira_ppanggolin.sh"

# Load Library 
library(tidyverse)

###########################################################
# PPanGGOLiN
ppanggolin_matrix <- read.delim(file = "matrix.csv", header = T, sep = ",", check.names = F)

###########################################################
# EggNOG-mapper
eggnog_mapper <- read.delim(file = "annotation/eggnogmapper/all_protein_families.emapper.annotations.tsv", header = T, sep = "\t", check.names = F)
colnames(eggnog_mapper)[colnames(eggnog_mapper) != "#query"] <- paste0("eggnog_", colnames(eggnog_mapper)[colnames(eggnog_mapper) != "#query"])
ppanggolin_matrix <- merge(x=ppanggolin_matrix, y=eggnog_mapper, by.x="Gene", by.y="#query", all.x=T, sort=F)

###########################################################
# COGclassifier
cogclassifier <- read.delim(file = "annotation/cogclassifier/cog_classify.tsv", header = T, sep = "\t", check.names = F)
colnames(cogclassifier)[colnames(cogclassifier) != "QUERY_ID"] <- paste0("cog_", colnames(cogclassifier)[colnames(cogclassifier) != "QUERY_ID"])
ppanggolin_matrix <- merge(x=ppanggolin_matrix, y=cogclassifier, by.x="Gene", by.y="QUERY_ID", all.x=T, sort=F)

###########################################################
# dbCAN
dbcan <- read.delim(file = "annotation/dbcan/overview.tsv", header = T, sep = "\t", check.names = F)
colnames(dbcan)[colnames(dbcan) != "Gene ID"] <- paste0("dbcan_", colnames(dbcan)[colnames(dbcan) != "Gene ID"])
ppanggolin_matrix <- merge(x=ppanggolin_matrix, y=dbcan, by.x="Gene", by.y="Gene ID", all.x=T, sort=F)

###########################################################
# AMRFinder
amrfinder <- read.delim(file = "annotation/amrfinder/amrfinder.tsv", header = T, sep = "\t", check.names = F)
colnames(amrfinder)[colnames(amrfinder) != "Protein id"] <- paste0("amr_", colnames(amrfinder)[colnames(amrfinder) != "Protein id"])
ppanggolin_matrix <- merge(x=ppanggolin_matrix, y=amrfinder, by.x="Gene", by.y="Protein id", all.x=T, sort=F)

###########################################################
# VFDB
vfdb <- read.delim(file = "annotation/vfdb/vfdb_header.tsv", header = T, sep = "\t", check.names = F)
colnames(vfdb)[colnames(vfdb) != "qseqid"] <- paste0("vfdb_", colnames(vfdb)[colnames(vfdb) != "qseqid"])
ppanggolin_matrix <- merge(x=ppanggolin_matrix, y=vfdb, by.x="Gene", by.y="qseqid", all.x=T, sort=F)

###########################################################
# RGI
rgi <- read.delim(file = "annotation/rgi/rgi.txt", header = T, sep = "\t", check.names = F)
colnames(rgi)[colnames(rgi) != "ORF_ID"] <- paste0("rgi_", colnames(rgi)[colnames(rgi) != "ORF_ID"])
ppanggolin_matrix <- merge(x=ppanggolin_matrix, y=rgi, by.x="Gene", by.y="ORF_ID", all.x=T, sort=F)

###########################################################
# Save table
ppanggolin_matrix[is.na(ppanggolin_matrix)] <- ""
ppanggolin_matrix[ppanggolin_matrix=="NA"] <- ""
write.table(ppanggolin_matrix, "matrix_annotated.tsv", row.names=FALSE, col.names=TRUE, quote=TRUE, sep="\t")



