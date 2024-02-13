# Immunopeptidome-genomic-coordintes
Script to convert immunopeptidome files generated with SOPRANO to genomic coordinates 
Programs required: R
File required: GRCh38_transcripts_gene_protein.csv

### Step 1: Add protein ID to protein coordinates ###

# Set the directory where your files are located
directory_path <- "/path/to/your/immunopeptidomes/directory"

# Get a list of files matching the pattern
file_list <- list.files(path = directory_path, pattern = "NAME1.bed", full.names = TRUE)

# Loop through each file
for (file_path in file_list) {
  # Read the immunopeptidome file
  immunopeptidome <- read.delim(file_path, header = FALSE)
  colnames(immunopeptidome) <- c("ENST", "Start", "End")
  
  # Read the reference file
  ref <- read.csv("GRCh38_transcripts_gene_protein.csv", sep = "\t", header = FALSE)
  colnames(ref) <- c("CHR", "Start", "End","ENSG", "ENST", "Symbol", "ENSP")
  
  # Fuse dataframes by the common column
  ENSP <- rep(NA, nrow(immunopeptidome))
  immunopeptidome <- data.frame(immunopeptidome, ENSP)
  immunopeptidome$ENSP <- ref$ENSP[match(immunopeptidome$ENST, ref$ENST)]
  
  # Write the modified data frame to a new file
  output_file_path <- paste0(sub(".bed$", "_proteinID.bed", basename(file_path)), sep = "")
  write.table(immunopeptidome, file = output_file_path, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("File processed:", file_path, "\n")
  cat("Output written to:", output_file_path, "\n\n")
}

### Step 2: conversion of protein coordinates to genomic coordinates ###

# Load libraries
library(ensembldb)
library(EnsDb.Hsapiens.v86)  #GRCh38
library(GenomicRanges)

# Load the annotation database containing the annotations from Ensembl release 86
edb <- EnsDb.Hsapiens.v86

# Define the pattern for file names
file_pattern <- "NAME2.bed"

# Get a list of files matching the pattern
files <- list.files(pattern = file_pattern)

# Create a loop to process multiple files
for (file_name in files) {
  # Read the data from the current file
  df <- read.delim(file = file_name, sep = "\t", header = TRUE)
  
  # Remove rows with NA values in the ENSP column
  df <- df[complete.cases(df$ENSP), ]
  
  # Check for missing or NA values in the ENSP column
  if (any(is.na(df$ENSP))) {
    stop(paste("Error: Missing or NA values found in the ENSP column of file", file_name))
  }
  
  # Define the IRanges object
  protein_cr <- IRanges(start = df$Start, end = df$End, names = df$ENSP)
  
  # Convert protein coordinates to genomic coordinates
  genomic_cr <- proteinToGenome(protein_cr, edb)
  
  # Define a function to filter ranges where the interval is greater than 27 nt to get rid of annotation that do not correspond to the length of the peptide
  filter_ranges <- function(x) { # this function takes the argument 'x' => GRanges object containing the genomic coordinates 
    x[width(x) > 27] # filter the range based on their width
  }
  
  # Apply the filter function to each protein's ranges
  genomic_cr_filtered <- lapply(genomic_cr, filter_ranges) # apply the 'filter_ranges" function to each element contained in the list 'genomic_cr'
  
  # Combine the filtered ranges into a single GRangesList
  genomic_cr_filtered <- unlist(genomic_cr_filtered) # take a list as an input and returns a single vector by concatenating the elements of the list

  # Create a data frame that contains the genomic coordinates, transcript_id, and the protein coordinates
  results_list <- lapply(seq_along(genomic_cr_filtered), function(j) {
    data.frame(
      seqnames = seqnames(genomic_cr_filtered[[j]]),
      start = start(genomic_cr_filtered[[j]]),
      end = end(genomic_cr_filtered[[j]]),
      tx_id = mcols(genomic_cr_filtered[[j]])$tx_id,
      protein_id = mcols(genomic_cr_filtered[[j]])$protein_id,
      protein_start = mcols(genomic_cr_filtered[[j]])$protein_start,
      protein_end = mcols(genomic_cr_filtered[[j]])$protein_end,
      stringsAsFactors = FALSE
    )
  })
  
  # Combine the list of data frames into a single data frame
  results <- do.call(rbind, results_list)
  colnames(results) <- c("Chr", "Chr_start", "Chr_end", "Tx_id", "Protein_id", "Protein_start", "Protein_end")
    
  # Save the results in a BED file
  output_file <- paste0(sub(".bed$", "_genomic_cr.bed", basename(file_name)), sep = "")
  write.table(results, file = output_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}
