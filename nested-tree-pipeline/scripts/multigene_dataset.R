#--------------------------------------------------------------------------------------
# Script Title: Single Gene Alignment and Multigene Construction
# Author: M.A. Thanuja M. Fernando
# Date: 2025-01-18
# Version: 1.0
# Description: 
# This script processes single gene data in fasta format (nucleotide sequence data) 
# Users are required to modify the .fas file name according to their specific dataset whenever new data (species or adding new genes) is available.
# This script aligns each gene individually. Performs multiple sequence alignment (MSA) 
# Finally, the aligned gene sequences are concatenated into a multi-gene sequence dataframe and converted into FASTA format for further analysis.
# Input Files: 4C4New.fas (nucleotide sequences in FASTA format) [change this file name to match your new dataset]
# Output Files: FinalAlignmentChap2_12193New.fasta (final concatenated multi-gene sequence dataset)  [you can modify this file name as per your preference]
# Contact: [thanuja@uoguelph.ca]


##This code is to clean the single genes in multiple gene dataset and combine single genes after adding new sequences and align each gene separately. Use these aligned genes and combine them into one dataset and convert it to fatsa file. 

#Load packages
library(seqinr)
library(DECIPHER)
library(Biostrings)
library(tidyr)
library(dplyr)
library(readr)
#read ncbi downloaded fasta in to R to see if there are new species in this list
fas_4c4_ncbi <-readDNAStringSet("4C4New.fas") #[change this file name to match your new single gene dataset which includes new species]

#convert the sequences into a data frame
sequences_df <- data.frame(head= names(fas_4c4_ncbi),sequence = as.character(fas_4c4_ncbi), stringsAsFactors = FALSE)
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
print(cleaned_data)

nchar(cleaned_data$sequence)

class(cleaned_data$head)
# Replace spaces with underscores for all entries in the 'head' column
cleaned_data$head <- gsub(" ", "_", cleaned_data$head)

# Find entries in 'head' column containing "_sp"
sp_entries <- cleaned_data$head[grep("_sp\\.$|_cf\\.$|_cf\\._|_aff\\.$|_sp\\._|_gen\\.$|_gr\\.$", cleaned_data$head)]
#check what are they
sp_entries
sp_entries_vec <- c(sp_entries)
#Now I check one by one useing Blast to get the epithets names of that genus.
#There's no specific name. they are unknown. Since, we need full species name. not genus, we remove these unknown species from our list. 
# Remove rows based on species names
cleaned_data <- cleaned_data[!cleaned_data$head %in% sp_entries_vec, ]
#to remove row names as a column.
rownames(cleaned_data) <- NULL
Final_4c4_alignment <- cleaned_data 
#Check details
nchar(cleaned_data$sequence)
class(cleaned_data)
dim(cleaned_data)
#Write a fasta file from final alignmnet
# Convert sequences to DNAStringSet object
sequences <- DNAStringSet(cleaned_data$sequence)

# Assign names to the sequences
names(sequences) <- cleaned_data$head

# Write the DNA sequences to a FASTA file
writeXStringSet(sequences, file = "Final4c4Alignment.fasta")
#-----------------------------------

# Path to the directory containing FASTA files
fasta_dir <- "/home/thanu/Desktop/FishData/Chap2/MultiGene/Thanu/Aligned" #[change this to your path]

# List all FASTA files in the directory
fasta_files <- list.files(path = fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

# Check if there are at least 27 files
if (length(fasta_files) < 27) {
  stop("You need at least 27 FASTA files.")
}

# Select the first 27 files
fasta_files <- fasta_files[1:27] #[If needed to add new genes you can modify this to number of genes you have]

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

#------------------------------
#First merge all 27 dataframes into a one datframe using seq.name
#multigene_4856 <- merge(fasta.12S,fasta.16S,fasta.4c4,fasta.COI,fasta.CytB,fasta.Enc1,fasta.Ficd,fasta.Glyt,fasta.Hoxc6a,fasta.Kiaa1239,fasta.Myh6,fasta.Nd2,fasta.ND4,fasta.Panx2,fasta.Plagl2,fasta.Ptr,fasta.Rag1,fasta.Rag2,fasta.Rhodopsin,fasta.Ripk4,fasta.Sh3px3,fasta.Sidkey,fasta.Sreb2,fasta.Svep1,fasta.Tbr,fasta.Vcpip,fasta.Zic1,by="seq.name")
# Join multiple data.frames
library(tidyverse)
#colnames( MCytB)
# Assuming you have already created and loaded your individual dataframes M12S, M16S, M4c4, etc.

# List of dataframes with adjusted names [add all your gene names]
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
filtered_df <- merged_df1[apply(merged_df1[,-1], 1, count_non_na) >= 2, ]

#Here's an extra step for counting the number of species in each gene.
#use the sapply function to count non-NA values in each gene column
#This will edxclude the "Species_name" column
gene_counts <- sapply(filtered_df [,-which(names(filtered_df)=="Species_name")],function(x)sum(!is.na(x)))
#> gene_counts
#M12S       M16S       M4c4       MCOI      MCytB      MEnc1      MFicd      MGlyt    MHoxc6a  MKiaa1239      MMyh6       MNd2       MND4     MPanx2 
#4619       6010       1313      10092       6631       1705        606       1055        122        851       2251       3531       2910        766 
#MPlagl2       MPtr      MRag1      MRag2 MRhodopsin     MRipk4    MSh3px3    MSidkey     MSreb2     MSvep1       MTbr     MVcpip      MZic1 
#1490       1945       4421       2129       3230        600       1455        926       1432        280       1815        374       2003 


#Extra step ends here
#Count number of markers per species
marker_count <- filtered_df %>% 
  dplyr::mutate(marker_count =rowSums(!is.na(filtered_df[,-1]))) %>% 
  dplyr::select(Species_name,marker_count)
write.csv(marker_count,"marker_count.csv")
min_markers <- min(marker_count$marker_count)
max_markers <- max(marker_count$marker_count)
avg_markers <- mean(marker_count$marker_count)

#Identify species with maximum markers (27 markers)
speciesWith27Markers <- marker_count %>% 
  filter(marker_count==27) %>% 
  pull(Species_name)

# Define a function to replace NA with dashes based on column size
replace_na_with_dashes <- function(column) {
  na_indices <- is.na(column)
  non_na_values <- column[!na_indices]
  
  if (length(non_na_values) > 0) {
    non_na_length <- nchar(non_na_values[1])
    dashes <- paste(rep("-", non_na_length), collapse = "")
    column[na_indices] <- dashes
  }
  
  return(column)
}

# Apply the function to each column in merged_df except ID
for (col in names(filtered_df)[-1]) { # Skip ID column
  filtered_df[[col]] <- replace_na_with_dashes(filtered_df[[col]])
}

nchar(filtered_df[5,])

# Check if each sequence in a column is aligned
check_alignment <- function(column) {
  lengths <- nchar(column)
  all(lengths == lengths[1])
}

# Apply the function to each column except ID and store results in a list
alignment_check <- sapply(filtered_df[-1], check_alignment)

# Combine columns from Sequence_12S to Sequence_COI into one column
filtered_df$Combined_Sequences <- apply(filtered_df[,-1], 1, paste, collapse = "")

# Keep only ID and Combined_Sequences columns
final_df <- filtered_df[, c("Species_name", "Combined_Sequences")]


# Check for duplicates in the 'Sequence' column
duplicates <- duplicated(final_df$Species_name) | duplicated(final_df$Species_name, fromLast = TRUE)

# Filter the dataframe to show rows where 'Sequence' column has duplicates
duplicate_rows <- final_df[duplicates, ]
# Keep only the first occurrence of each unique sequence
df_unique <- final_df[!duplicated(final_df$Species_name), ]
# Check for duplicates in the 'Sequence' column
duplicates <- duplicated(df_unique$Combined_Sequences) | duplicated(df_unique$Combined_Sequences, fromLast = TRUE)

# Filter the dataframe to show rows where 'Sequence' column has duplicates
duplicate_rows <- df_unique[duplicates, ]
# Keep only the first occurrence of each unique sequence
df_uniqueNew <- df_unique[!duplicated(df_unique$Combined_Sequences), ]
nchar(final_df$Combined_Sequences)

#Function to count nucleotide based for each species
count_bases <- function(sequence){
  nchar(gsub("-","",sequence))
}
#Appply the function to cont nucleotide based for each species
NucleotideBAsesCount <- df_uniqueNew %>% 
  mutate(base_count=sapply(Combined_Sequences,count_bases))

#Calculate the min,max and avg number of nucleotides per species (without gaps)
min_bases <- min(NucleotideBAsesCount$base_count)
max_bases <- max(NucleotideBAsesCount$base_count)
avarage_bases <-mean(NucleotideBAsesCount$base_count)  

#Write a fasta file from final alignmnet
# Convert sequences to DNAStringSet object
sequences <- DNAStringSet(df_uniqueNew$Combined_Sequences)

# Assign names to the sequences
names(sequences) <- df_uniqueNew$Species_name

# Write the DNA sequences to a FASTA file
writeXStringSet(sequences, file = "FinalAlignmentChap2_12193New.fasta") #[you can modify this file name as per your preference]
