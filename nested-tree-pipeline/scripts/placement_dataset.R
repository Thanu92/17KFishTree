#This code is to prepare the COI abrcoding sequences dataset in fasta format to place in the multigene backbone tree to have the placement tree
#This should run with multigene_dataset.R 
#------------------------------------------------------------------
# Script Title: COI Barcoding Alignment and Placement Preparation
# Author: M.A. Thanuja M. Fernando
# Date: 2025-01-18
# Version: 1.1
# Description: 
# This script processes COI barcoding sequences in FASTA format to prepare them for placement in a multigene backbone tree.
# Users should update the .fas file name according to their dataset when adding new species or genes.
# The script aligns COI sequences, removes incomplete or unidentified species, and ensures data consistency for phylogenetic placement.
# Input Files: COI.Mask.fst (COI barcoding sequences in FASTA format) [update with new COI sequences]
# Output Files: df_COI_5218_new.fasta (final aligned COI dataset)  [modify as needed]
# Contact: [thanuja@uoguelph.ca]
# This script should be run alongside multigene_dataset.R

#Load packages
library(Biostrings)
library(tidyverse)

#read ncbi downloaded fasta in to R to see if there are new species in this list
fas_coi_ncbi <-readDNAStringSet("COI.Mask.fst") ##[update this file name with COI sequences for new species]

#convert the sequences into a data frame
sequences_df <- data.frame(head= names(fas_coi_ncbi),sequence = as.character(fas_coi_ncbi), stringsAsFactors = FALSE)
#Need to extract species name
#remove the word "PREDICTED" from head column
sequences_df$head <- gsub("PREDICTED: ", "",sequences_df$head)

# Function to extract species name

extract_species_name <- function(head) {
  match <- regexpr("\\b[A-Z][a-z]+ [a-z]+\\b", head)
  if (match[1] != -1) {
    species_name <- regmatches(head, match)
    return(species_name)
  } else {
    return(head)  # If no match is found, return the original description
  }
}

# Use mutate to apply the function to the description column
cleaned_data <- sequences_df %>%
  mutate(head = sapply(head, extract_species_name))

# Print the result
#print(cleaned_data)

nchar(cleaned_data$sequence)

class(cleaned_data$head)
# Replace spaces with underscores for all entries in the 'head' column
cleaned_data$head <- gsub(" ", "_", cleaned_data$head)

# Find entries in 'head' column containing "_sp"
sp_entries <- cleaned_data$head[grep("_sp\\.$|_cf\\.$|_cf\\._|_aff\\.$|_sp\\._|_gen\\.$|_gr\\.$|_aff\\._", cleaned_data$head)]
#check what are they
sp_entries
sp_entries_vec <- c(sp_entries)
#Now I check one by one useing Blast to get the epithets names of that genus.
#There's no specific name. they are unknown. Since, we need full species name. not genus, we remove these unknown species from our list. 
# Remove rows based on species names
cleaned_data <- cleaned_data[!cleaned_data$head %in% sp_entries_vec, ]
#to remove row names as a column.
rownames(cleaned_data) <- NULL
Final_CO1_alignment <- cleaned_data 
#Check details
nchar(cleaned_data$sequence)
class(cleaned_data)
dim(cleaned_data)
names(cleaned_data)
df_unique <- cleaned_data %>%
  distinct(head, .keep_all = TRUE)

cleaned_data <- df_unique
#Write a fasta file from final alignmnet
# Convert sequences to DNAStringSet object
sequences <- DNAStringSet(cleaned_data$sequence)

# Assign names to the sequences
names(sequences) <- cleaned_data$head

# Write the DNA sequences to a FASTA file
writeXStringSet(sequences, file = "FinalCO1alignment.fasta")
#-----------------------------------------------------
#Check coi sequences not in multi-gene dataset to prepare for as the placement sequences
# Path to the directory containing FASTA files
fasta_dir <- "/home/thanu/Desktop/FishData/Chap2/MultiGene/Thanu/Aligned" #[change this to your path]

# List all FASTA files in the directory
fasta_files <- list.files(path = fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

# Check if there are at least 27 files
if (length(fasta_files) < 27) {
  stop("You need at least 27 FASTA files.")
}

# Select the first 27 files
fasta_files <- fasta_files[1:27]

# Function to convert a FASTA file to a dataframe and remove row names as a column
fasta_to_dataframe <- function(file) {
  fasta <- readDNAStringSet(file)
  df <- data.frame(
    ID = names(fasta),
    Sequence = as.character(fasta),
    stringsAsFactors = FALSE
  )
  # Remove row names as a column
  rownames(df) <- NULL
  return(df)
}

# Extract base names without extension
file_names <- tools::file_path_sans_ext(basename(fasta_files))

# Loop through each file and create a dataframe, assigning it to the global environment
for (i in 1:length(fasta_files)) {
  df_name <- paste0("M", file_names[i])  # Prefix with "M"
  assign(df_name, fasta_to_dataframe(fasta_files[i]), envir = .GlobalEnv)
}

#----
#First merge all 27 dataframes into a one datframe using seq.name
# you have already created and loaded your individual dataframes M12S, M16S, Mcoi, etc.

# List of dataframes with adjusted names
list_df <- list(M12S, M16S, M4c4, MCOI, MCytB, MEnc1, MFicd, MGlyt, MHoxc6a, 
                MKiaa1239, MMyh6, MNd2, MNd4, MPanx2, MPlagl2, MPtr, 
                MRag1, MRag2, MRhodopsin, MRipk4, MSh3px3, MSidkey, MSreb2, 
                MSvep1, MTbr, MVcpip, MZic1)

# Define a function to merge two data frames by ID
merge_by_id <- function(df1, df2) {
  merge(df1, df2, by = "ID", all = TRUE)
}

# Use Reduce to merge all data frames in the list
merged_df <- Reduce(merge_by_id, list_df)

#Remove first three rows
# Rearrange columns and rename if necessary
merged_df1 <- merged_df[-(1:3),]
nchar(merged_df1[5,])

# rename
colnames(merged_df1) <- c("Species_name", "M12S", "M16S", "M4c4", "MCOI", "MCytB", "MEnc1", "MFicd", "MGlyt", "MHoxc6a","MKiaa1239", "MMyh6", "MNd2", "MND4", "MPanx2","MPlagl2", "MPtr", "MRag1", "MRag2", "MRhodopsin", "MRipk4", "MSh3px3", "MSidkey", "MSreb2", "MSvep1", "MTbr", "MVcpip", "MZic1")
merged_df1_wtNames <- merged_df1


# Filter rows to keep only those with two or more non-NA columns
# Count non-NA values in each row
count_non_na <- function(row) {
  sum(!is.na(row))
}

# # Apply the function to filter rows
filtered_df1 <- merged_df1[apply(merged_df1[,-1], 1, count_non_na) >= 2 & !is.na(merged_df1$MCOI), ]
# # Apply the function to filter rows
filtered_df <- merged_df1[apply(merged_df1[,-1], 1, count_non_na) >= 2, ]
names(filtered_df1)
difference <- anti_join(filtered_df, filtered_df1, by = "Species_name")

not_in_co1 <-c(difference$Species_name) 

#Check the difference between MCO1  15403 and filtered df...
#col names should be the same
MCOI_new <- MCOI
colnames(MCOI_new)[1] <- "Species_name"
COI_difference <- anti_join(MCOI_new, filtered_df, by = "Species_name") %>%
  select(names(MCOI_new))

nchar(MCOI_new$Sequence)

# Check for duplicates in the 'Sequence' column
duplicates <- duplicated(COI_difference$Species_name) | duplicated(COI_difference$Species_name, fromLast = TRUE)

# Filter the dataframe to show rows where 'Sequence' column has duplicates
duplicate_rows <- COI_difference[duplicates, ] #NO duplicates in species names for COI
# # Keep only the first occurrence of each unique sequence
# df_unique <- final_df[!duplicated(final_df$Species_name), ]
# Check for duplicates in the 'Sequence' column
duplicates <- duplicated(COI_difference$Sequence) | duplicated(COI_difference$Sequence, fromLast = TRUE)

# Filter the dataframe to show rows where 'Sequence' column has duplicates
duplicate_rows <- COI_difference[duplicates, ]
# Keep only the first occurrence of each unique sequence
df_COI_barcode <- COI_difference[!duplicated(COI_difference$Sequence), ]
nchar(COI_difference$Sequence)

# Check if each sequence in a column is aligned
check_alignment <- function(column) {
  lengths <- nchar(column)
  all(lengths == lengths[1])
}

# Apply the function to each column except ID and store results in a list
alignment_check <- sapply(df_COI_barcode[-1], check_alignment)


#Function to count nucleotide based for each species
count_bases <- function(sequence){
  nchar(gsub("-","",sequence))
}
names(df_COI_barcode)
#Appply the function to cont nucleotide based for each species
NucleotideBAsesCount <- df_COI_barcode %>% 
  mutate(base_count=sapply(Sequence,count_bases))

#Calculate the min,max and avg number of nucleotides per species (without gaps)
min_bases <- min(NucleotideBAsesCount$base_count)
max_bases <- max(NucleotideBAsesCount$base_count)
avarage_bases <-mean(NucleotideBAsesCount$base_count)

#Write a fasta file from final alignmnet
# Convert sequences to DNAStringSet object
sequences <- DNAStringSet(df_COI_barcode$Sequence)

# Assign names to the sequences
names(sequences) <- df_COI_barcode$Species_name

# Write the DNA sequences to a FASTA file
writeXStringSet(sequences, file = "df_COI_5218_new.fasta") #[you can modify this file name as per your preference]
