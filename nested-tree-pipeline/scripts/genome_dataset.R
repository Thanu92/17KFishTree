#----------------------------------------------------------------------
# Script Title: BUSCO Gene Alignment and Genome Construction from NCBI .faa Files
# Author: M.A. Thanuja M. Fernando
# Date: 2025-01-12
# Version: 1.0
# Description: 
# This script processes BUSCO gene data downloaded from NCBI in .faa format (protein sequence data), where all genes and species are listed together. 
# Users are required to modify the .faa file name according to their specific dataset whenever new data is available.
# The script first organizes and sorts the data, ensuring each gene's sequences are correctly associated with their respective species. 
# After sorting, the script aligns each BUSCO gene individually. Performs multiple sequence alignment (MSA) 
# Finally, the aligned gene sequences are concatenated into a single genome sequence and converted into FASTA format for further analysis.
# It handles protein FASTA files, generates individual gene files, and writes out the final aligned genome.
# Input Files: protein.faa (protein sequences in FASTA format) [change this file name to match your new dataset]
# Output Files: genomeNew.fasta (final concatenated genome sequence)  [you can modify this file name as per your preference]
# Contact: [thanuja@uoguelph.ca]

#The code to generate final fasta file from aligned protein sequences of busco genes to prepare genome fasta data file.
#Load packages
require(tidyverse)
require(tidyr)
library(ampir)#to read .faa files in r
library(dplyr)
library(Biostrings)
library(foreach)
# if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
# BiocManager::install("msa")
library(msa)
library(stringr)
getwd()
#Read protein fasta files to r downloaded from ncbi
pro_fasta <- read_faa(file="protein.faa")

gene_filter <- function(protein_set){
  new <- separate(protein_set, col = seq_name , into = paste0("value",1:6),sep=" ")
  new1 <- new[,c(2,3,4,6)]
  new1 <- cbind(new1,protein_set$seq_aa)
  
  #Combine value 3 and 4 together (genus and species anme in two different columns)
  new2 <- new1 %>%
    unite(species_name, c("value3", "value4"))
  # Extract strings after "="
  new2$species_name <- str_extract(new2$species_name, "(?<==).*")
  #remove the back bracket in species column
  # Remove the symbol from dataframe values
  new2$species_name <- str_replace(new2$species_name, "]", "")
  
  #Do the same thing to remove [] in value6 column
  
  # Extract strings after "="
  new2$value6 <- str_extract(new2$value6, "(?<==).*")
  #remove the back bracket in species column
  # Remove the symbol from dataframe values
  new2$value6 <- str_replace(new2$value6, "]", "")
  
  #rename columns
  colnames(new2) <- c("gene","species_name","isoform","seq_aa")
  
  #When adding outgroup thene gene names are in uppercase. R is case sensitive. 
  #Hence, same gene name in upper and lower case act as two genes. So, lets make all in lower case
  new2$gene <- tolower(new2$gene)
  
  #eaxtract the fiest two columns and combine first two columns in to one, then we can easily remove duplicated rows
  new_Df <- new2[,c(1:2)]
  new_Df$Colnew <- apply(new_Df, 1, function(x) paste(na.omit(x), collapse = " "))
  
  #Combine the colnew to the new2 dataframe, we need the sequences as well
  new_df <- cbind(new_Df$Colnew,new2)
  names(new_df)
  
  #remove isofrom X2:4 and keep only isoform X1 sequences
  df1 <-
    new_df[!grepl("X2|X3|X4|X5|X6|X7|X8|X9", new_df$isoform),]
  #Get unique rows. we need only isoform x1
  unique_df <- df1[!duplicated(df1$`new_Df$Colnew`), ]
  
  #get a subdataframe
  gene_speciec_Df <- unique_df[,c(2,3,5)]
  
  #remove rows with gene with LOC^ (coz it's not a gene name)
  gene_speciec_Df_new <- gene_speciec_Df[!grepl('^LOC|^lg|orf|^si$|^si:',gene_speciec_Df$gene),]
  
  #divide dataframes into several dataframe using gene name
  gene <- gene_speciec_Df_new %>% group_by(gene)
  
  gene_split <- group_split(gene)
  gene_split[[8]]
  #Remove dataframes in the list that have less than 3 including 3 (just in case if these are only from three out groups)
  gene_split <- gene_split[sapply(gene_split, nrow) > 3]
  return(gene_split)
  
}

#### ----------------------------------STEP2

#This step is to get individual dataframes from the list of dataframes

#Function for msa
MSA_genes <- function(gene_split){
  #This step is to get individual dataframes from the list of dataframes
  
  # Commenting this part out for now
  # list2env(setNames(gene_split,paste0("newdf",seq_along(gene_split))),envir = .GlobalEnv)
  
  # Instead we could do this
  
  # first for each df, I will concatenate gene with species name using paste for the header info foreach fasta file
  gene_split <- foreach(i=1:length(gene_split)) %do% {
    gene_split[[i]] %>% mutate(header = paste0(gene_split[[i]]$gene, ";", gene_split[[i]]$species_name, sep=""))
  }
  
  # Then we can just create a new directory for the fasta files
  dir <- "fastaFiles"
  dir.create(dir, showWarnings = FALSE)
  
  # Convert each dataframe to aaStringSet format
  aaStringSetList <- foreach(i=1:length(gene_split)) %do% { AAStringSet(x = as.character(gene_split[[i]]$seq_aa)) }
  
  # Name each AAStringSet
  foreach(i=1:length(gene_split)) %do% { aaStringSetList[[i]]@ranges@NAMES <- as.character(gene_split[[i]]$header) }
  
  # Create filenames for each fasta, can name by gene name
  filenames <- foreach(i=1:length(gene_split)) %do% { paste0(dir, "/", gene_split[[i]]$gene[1], ".fas", sep="") }
  
  # Then use writeXStringSet to write out each one to a fasta file in the fastaFiles directory
  foreach(i=1:length(gene_split)) %do% writeXStringSet(aaStringSetList[[i]], filenames[[i]], format = "fasta")
  
  # Then to align each one you could do (just using msa package for example)
  
  
  # Grab each fasta file and import back in
  fastaFiles <- foreach(i=1:length(gene_split)) %do% readAAStringSet(filenames[[i]], format = "fasta")
  
  # Perform an msa on each file (may take a while to complete each one)
  # to test it works: alignments <- foreach(i=1:2) %do% msa(fastaFiles[[i]], verbose=TRUE)
  alignments <- foreach(i=1:length(gene_split)) %do% msa(fastaFiles[[i]], verbose=FALSE) # If you want to see the progress for every alignment, set to TRUE
  alignments[[1]]
  #alignments_100_1 <- foreach(i=1:length(gene_split)) %do% msa(fastaFiles[[i]],method =c("Muscle"),type="protein", verbose=FALSE) # If you want to see the progress for every alignment, set to TRUE
  class(alignments[[1]])
  
  # Then we can just create a new directory for the fasta files
  dir <- "AlnfastaFiles"
  dir.create(dir, showWarnings = FALSE)
  # Convert each dataframe to aaStringSet format
  AlnaaStringSetList <- foreach(i=1:length(alignments)) %do% { AAStringSet(x = as.character(alignments[[i]])) }
  
  # Create filenames for each fasta, can name by gene name
  Alnfilenames <- foreach(i=1:length(gene_split)) %do% { paste0(dir, "/", gene_split[[i]]$gene[1], ".fas", sep="") }
  
  # Then use writeXStringSet to write out each one to a fasta file in the fastaFiles directory
  foreach(i=1:length(alignments)) %do% writeXStringSet(AlnaaStringSetList[[i]], Alnfilenames[[i]], format = "fasta")
  
  
}
protein_set <- gene_filter(pro_fasta) #[change this to match your new fasta file]
MSA_genes(protein_set)

#--------------------------

#Now we need to concatanate all the genes together
#First read all the fas fies into r

files <- list.files(path="C:/Users/Thanu Fdo/OneDrive/Desktop/FishData_2022/Chap2/Genome/AlnfastaFiles/",pattern="*.fas",full.names=T,recursive=F)
fastaList <- foreach(i=1:length(files))%do%read_faa(files[i])
first <- fastaList[[1]]
#Convert to a list of dataframes

fastaList_new <- foreach(i=1:length(fastaList))%do%as.data.frame(fastaList[[i]])

#split the column seq_name using ';' as the separator in to gene and species names.Since, all the dataframes have the same column use lapply,
fastaList_species_names <- lapply(fastaList_new, tidyr::separate, col = "seq_name", into = c("gene", "species_name"), sep = ";")
names(fastaList_species_names[[1]])#"gene" "species_name" "seq_aa" 
check <- fastaList_species_names[[2]]
class(fastaList_species_names)
#Since we need the species and aa sequences we get the 2nd and 3rd column 
fastaList_species_names_new <- lapply(fastaList_species_names , FUN = function(x){x[, c(2:3)]})
names(fastaList_species_names_new[[1]])
#Horizontally concatanate the multiple dataframes into a one dataframe using species_name column
big_df <- Reduce(function(x, y) merge(x, y, by = "species_name", all = TRUE), fastaList_species_names_new)
#warnings(big_df)

# #get the gene names
#Since we need the gene column we get the 1st column of fastaList_species_names dataframe
gene_names <- lapply(fastaList_species_names , FUN = function(x){x[,1]})

# Using lapply to extract the first element from each data frame in the list
first_elements <- lapply(gene_names, function(x) x[[1]])
#get the list elements as a vector
vector_elements <- unlist(first_elements)
#rename gene names in the dataframe "genome_df" using the vector_elements vector

# Dynamically assign column names based on vector_elements
colnames(big_df)[2:length(big_df)] <- vector_elements

# Replace NAs with dashes
replace_na_with_dashes <- function(x) {
  non_na_values <- x[!is.na(x)]
  if (length(non_na_values) > 0) {
    dashes <- paste(rep("-", nchar(non_na_values[1])), collapse = "")
    x[is.na(x)] <- dashes
  }
  return(x)
}

cols_to_process <- vector_elements
big_df[cols_to_process] <- lapply(big_df[cols_to_process], function(x) replace_na_with_dashes(x))

# Merge all genes dynamically
first_gene <- vector_elements[1]
last_gene <- vector_elements[length(vector_elements)]
NewGenome_dfSep <- big_df %>%
  unite("Merged", first_gene:last_gene, sep="")

# Write the final genome to FASTA
df_to_faa(NewGenome_dfSep, file = "genomeNewtest.fasta") #[you can modify this file name as per your preference]
