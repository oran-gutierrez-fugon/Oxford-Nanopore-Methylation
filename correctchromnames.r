library(data.table)
library(foreach)
library(doParallel)

# Set the number of cores/workers to utilize
num_cores <- 20  # Adjust this based on your system's capabilities

# Initialize parallel backend
registerDoParallel(cores = num_cores)  # Use doParallel, or choose appropriate backend

# Read the dictionary file
dictionary_file <- "/share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/hg38refGenome/GCF_000001405.39_GRCh38.p13_assembly_report.chrnames"
dictionary <- fread(dictionary_file, header = TRUE, colClasses = c("character", "character"))

# Create a named vector from the dictionary
chromosome_mapping <- setNames(dictionary$DesiredName, dictionary$CurrentName)

# Path to your BED or BEDGRAPH file
bed_file <- "/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/cpg/NP4-4cpg.bed"

# Path to output the modified BED/BEDGRAPH file
output_file <- "/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/cpg/NP4-4cpgMODIFIED.bed"


# Read the BED/BEDGRAPH file in chunks
chunk_size <- 100000  # Adjust the chunk size as per your system's capabilities
reader <- data.table::fread(bed_file, sep = "\t", header = FALSE, colClasses = "character")

# Define a function to modify chromosome names
modify_chromosome_names <- function(chunk) {
  chunk$V1 <- chromosome_mapping[chunk$V1]
  return(chunk)
}

# Apply parallel processing to modify chromosome names
modified_data <- foreach(chunk = iter(reader, by = chunk_size), .combine = rbind) %dopar% {
  modify_chromosome_names(chunk)
}

# Write the modified BED/BEDGRAPH file
fwrite(modified_data, output_file, sep = "\t", quote = FALSE)

# Stop parallel backend
stopImplicitCluster()








library(data.table)
library(foreach)
library(doParallel)

# Set the number of cores/workers to utilize
num_cores <- 50  # Adjust this based on your system's capabilities

# Initialize parallel backend
registerDoParallel(cores = num_cores)  # Use doParallel, or choose appropriate backend

# Read the dictionary file
dictionary_file <- "/share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/hg38refGenome/GCF_000001405.39_GRCh38.p13_assembly_report.chrnames"
dictionary <- fread(dictionary_file, header = TRUE, colClasses = c("character", "character"))

# Create a named vector from the dictionary
chromosome_mapping <- setNames(dictionary$DesiredName, dictionary$CurrentName)

# Path to your BED or BEDGRAPH file
bed_file <- "/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/cpg/NP4-4cpg.bed"

# Path to output the modified BED/BEDGRAPH file
output_file <- "/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/cpg/NP4-4cpgMODIFIED.bed"

# Read the BED/BEDGRAPH file in chunks
chunk_size <- 100000  # Adjust the chunk size as per your system's capabilities
reader <- data.table::fread(bed_file, sep = "\t", header = FALSE, colClasses = "character")

# Define a function to modify chromosome names
modify_chromosome_names <- function(chunk) {
  chunk$V1 <- chromosome_mapping[chunk$V1]
  return(chunk)
}

# Apply parallel processing to modify chromosome names
modified_data <- foreach(chunk = chunkIndex(reader, chunk_size), .combine = rbind) %dopar% {
  chunk_data <- reader[chunk]
  modify_chromosome_names(chunk_data)
}

# Stop parallel backend
stopImplicitCluster()

# Combine the modified chunks into a single data.table
modified_data <- rbindlist(modified_data)

# Write the modified BED/BEDGRAPH file
fwrite(modified_data, output_file, sep = "\t", quote = FALSE)




library(data.table)
library(foreach)
library(doParallel)

# Set the number of cores/workers to utilize
num_cores <- 30 # Adjust this based on your system's capabilities

# Initialize parallel backend
registerDoParallel(cores = num_cores)  # Use doParallel, or choose appropriate backend

# Read the dictionary file
dictionary_file <- "/share/lasallelab/Oran/dovetail/luhmes/methylation/fastqconcats/hg38refGenome/GCF_000001405.39_GRCh38.p13_assembly_report.chrnames"
dictionary <- fread(dictionary_file, header = TRUE, colClasses = c("character", "character"))

# Create a named vector from the dictionary
chromosome_mapping <- setNames(dictionary$DesiredName, dictionary$CurrentName)

# Path to your BED or BEDGRAPH file
bed_file <- "/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/cpg/NP4-4cpg.bed"

# Path to output the modified BED/BEDGRAPH file
output_file <- "/share/lasallelab/Oran/dovetail/luhmes/methylation/concatenated/NP4-4/cpg/NP4-4cpgMODIFIED.bed"

# Read the BED/BEDGRAPH file in chunks
chunk_size <- 100000  # Adjust the chunk size as per your system's capabilities
reader <- data.table::fread(bed_file, sep = "\t", header = FALSE, colClasses = "character")

# Define a function to modify chromosome names
modify_chromosome_names <- function(chunk) {
  chunk[, V1 := chromosome_mapping[match(V1, names(chromosome_mapping))]]
  return(chunk)
}

# Apply parallel processing to modify chromosome names
modified_data <- foreach(chunk = iter(reader, by = chunk_size), .combine = rbind) %dopar% {
  modify_chromosome_names(chunk)
}

# Stop parallel backend
stopImplicitCluster()

# Write the modified BED/BEDGRAPH file
fwrite(modified_data, file = output_file, sep = "\t", quote = FALSE)
