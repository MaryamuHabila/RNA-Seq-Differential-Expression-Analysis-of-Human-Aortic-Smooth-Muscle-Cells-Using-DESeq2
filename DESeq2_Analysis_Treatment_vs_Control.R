#Install packages
#install.packages("BiocManager")
#BiocManager::install("tximport")
#install.packages("readr")
#BiocManager::install("DESeq2")
#install.packages("tidyverse")
#install.packages("data.table")
#BiocManager::install("apeglm")

#Load libraries
library(tximport)
library(readr)
library(DESeq2)
library(data.table)
library(rtracklayer)


#Samples are grouped based on experimental conditions.(Creating a Sample Metadata Data Frame)
#A data.frame called samples is created and this stores metadata for the RNA-seq experiments. RNA-seq experiment was performed in Human aortic smooth muscle cells (Controls vs Treatment(cells treated with siRNA targeting the gene of interest)) Each row represents a sample, with two the columns:
#column1: sample_id: A unique identifier for each sample.
#colum 2: condition: A factor variable indicating whether the sample belongs to the "Control" or "Treatment" group.
samples <- data.frame(
  sample_id = c("1_21_Ctrl", "3_K1_Ctrl", "4_R1_Ctrl", "5_R1_Ctrl",
                "2_2SLT", "3_1SLT", "4_SLK1", "5_SLK2"),
  condition = factor(c(rep("Control", 4), rep("Treatment", 4)))
)


# Set row names as sample IDs 
#The row names of the samples data frame are set to the values in the sample_id column. This allows easier indexing and referencing of samples by their unique identifiers. After this step, the sample_id column is still present, but the row names provide a more convenient way to access the data
rownames(samples) <- samples$sample_id


#Define paths to quant.sf files (Salmon Quantification Files)
# Replace `/path/to/quant/files` with the directory containing the Salmon output
dir <- "/Users/maryamu/Documents/quantification_data_SMCs/quantification_data_SMCs"


# Constructing File Paths for Salmon Quantification Files: File paths are constructed using sample IDs and the new .quant.sf naming convention
#This step generates the full file paths to the quant.sf files for each sample by: Taking the base directory path stored in dir. Appending each sample_id from the samples data frame with the .quant.sf extension. Using file.path() to ensure correct path formatting across different operating systems. The resulting files vector contains the full paths to the quantification files, making it easier to read them into R for downstream analysis.
files <- file.path(dir, paste0(samples$sample_id, ".quant.sf"))


# Assigning Sample Names to File Paths 
#This step sets the names of the files vector to the corresponding sample_id values from the samples data frame. This makes it easier to reference specific samples when working with the files vector, as each file path is now labeled with its respective sample ID.For example, instead of accessing files by numeric index, you can retrieve a specific sample’s file path using: files["1_21_Ctrl"]
names(files) <- samples$sample_id


# Check that the file paths were created correctly
print(files)


# Verify that all files exist
file.exists(files)


#Read tx2gene Table The tx2gene table is used to map transcript IDs (e.g., those from the quant.sf files) to their corresponding gene IDs. This mapping is essential to aggregate transcript-level expression data into gene-level data

# Defining the Path to the GTF Annotation File
#This assigns the file path of a UCSC GTF (Gene Transfer Format) annotation file to the variable gtf_file. The GTF file contains gene annotations for the hg38 human genome, which are essential for transcript quantification, and feature counting.
gtf_file <- "/Users/maryamu/Desktop/hg38.refGene.gtf"


# Importing and Previewing the GTF Annotation File
gtf_data <- import(gtf_file) #This loads the GTF file using the rtracklayer package.The resulting object is a data.frame
head(gtf_data) #Shows the first few entries to verify successful import and inspect the structure of the annotation data.


#Creating a Transcript-to-Gene Mapping Data Frame from GTF Data
#This step extracts relevant genomic features from the imported GTF file (gtf_data) and stores them in a structured data.frame called tx2gene. The columns include: 1). Genomic coordinates: seqnames, ranges_start, ranges_end, strand; 2). Metadata from GTF attributes: source, type, score, phase; 3). Gene and transcript information: gene_id, transcript_id, gene_name. 4). Exon details: exon_number, exon_id
#This structured format is useful for RNA-seq analysis, particularly for transcript-level quantification and annotation.
tx2gene <- data.frame(
  seqnames = as.character(seqnames(gtf_data)),  # Extract seqnames
  ranges_start = start(ranges(gtf_data)),      # Extract start positions of ranges
  ranges_end = end(ranges(gtf_data)),          # Extract end positions of ranges
  strand = as.character(strand(gtf_data)),     # Extract strand
  source = mcols(gtf_data)$source,             # Metadata column
  type = mcols(gtf_data)$type,                 # Metadata column
  score = mcols(gtf_data)$score,               # Metadata column
  phase = mcols(gtf_data)$phase,               # Metadata column
  gene_id = mcols(gtf_data)$gene_id,           # Metadata column
  transcript_id = mcols(gtf_data)$transcript_id,  # Metadata column
  gene_name = mcols(gtf_data)$gene_name,       # Metadata column
  exon_number = mcols(gtf_data)$exon_number,   # Metadata column
  exon_id = mcols(gtf_data)$exon_id            # Metadata column
)



#Creating a Unique Gene Metadata Data Frame (Ensuring tx2gene has gene-level information)
#This step extracts unique gene entries from the tx2gene data frame and sets gene IDs as row names:
#The resulting gene_metadata table provides a concise mapping of genes with their associated annotations, useful for downstream analysis.
gene_metadata <- tx2gene[!duplicated(tx2gene$gene_id), ] #Removes duplicate gene_id entries, keeping only the first occurrence, creating a non-redundant gene_metadata table.
rownames(gene_metadata) <- gene_metadata$gene_id #Sets gene_id as row names for easy lookup.
head(tx2gene) #Displays the first few rows of tx2gene to inspect its structure.


#Import Quantification Data Using tximport
#The tximport step converts transcript-level data into gene-level data.#This step imports transcript-level quantifications from Salmon and aggregates them to the gene level using the tx2gene table.


#Extracting Gene Names from Quantification Files
#This code reads the first Salmon quantification file (files[1]) and extracts the Name column, which typically contains the gene identifiers. The extracted gene names are stored in the quant_ids vector.
quant_ids <- read.delim(files[1], header = TRUE)$Name #Reads the first quantification file as a tab-delimited file and selects the Name column (usually representing gene names or IDs).
head(quant_ids) #Displays the first few entries of the quant_ids vector to verify the data.



# Previewing Transcript IDs in the tx2gene Data Fram
#This step displays the first few transcript_id values from the tx2gene data frame. The transcript_id column contains identifiers for the transcripts associated with each gene.
head(tx2gene$transcript_id) #Shows the first few transcript IDs to verify the data and ensure that the correct transcript identifiers are included in the tx2gene data frame for further analysis.


# Finding Common Transcript IDs Between Quantification and Gene Dat
#This step finds the intersection of transcript IDs between the quant_ids vector (extracted from the quantification file) and the transcript_id column in the tx2gene data frame.This step ensures that the transcript IDs from the quantification files align with the ones in the gene metadata, enabling correct downstream analysis.
intersect(quant_ids, tx2gene$transcript_id) #Returns the transcript IDs that are present in both quant_ids (the list of IDs from the quantification data) and tx2gene$transcript_id (the list of IDs from the gene annotation data).


#Creating a Simplified Transcript-to-Gene Mapping Data Frame
#This step creates a tx2gene data frame that links transcript IDs to gene IDs. It simplifies the transcript_id by removing version numbers (e.g., "ENST00000123456.1" becomes "ENST00000123456") using the sub() function.This simplified tx2gene table helps with matching transcripts to genes during RNA-seq analysis
tx2gene <- data.frame(
  transcript_id = sub("\\..*$", "", mcols(gtf_data)$transcript_id),  #Strips version numbers from the transcript_id column in the gtf_data metadata.
  gene_id = mcols(gtf_data)$gene_id #Extracts the corresponding gene_id for each transcript.
)

#Extracting Transcript IDs from GTF Data
#This step retrieves the transcript_id column from the metadata (mcols()) of the gtf_data object. The transcript_id represents unique identifiers for each transcript in the GTF annotation file, which is crucial for mapping and analyzing transcript-level data.
mcols(gtf_data)$transcript_id #Extracts the transcript_id attribute, allowing you to associate each transcript with its corresponding gene and other features.


#Importing Transcript-Level Quantification Data with tximport
#This step uses the tximport() function to import transcript-level quantification data from Salmon (files) into R. The tx2gene mapping is provided to associate transcripts with genes, and the ignoreTxVersion = TRUE argument ensures version numbers are ignored during transcript matching.
#type = "salmon": Specifies the input data is from Salmon.
#tx2gene = tx2gene: Links transcript IDs to gene IDs using the previously defined tx2gene data frame.
#ignoreTxVersion = TRUE: Removes version numbers from transcript IDs to avoid mismatches.This prepares the data for downstream analysis, such as gene-level expression quantification.
txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene, 
                ignoreTxVersion = TRUE)


###### Run DESeq2 Analysis ######
######.   Comparison ######: 3 CONTROLS VS 3 TREATMENT GROUPS
#Subsetting Samples for Specific Conditions
#This step creates a new data frame samples_subset4 by selecting a subset of samples from the samples data frame. The samples chosen are identified by their sample_id values.
samples_subset4 <- samples[c("4_R1_Ctrl", "5_R1_Ctrl","1_21_Ctrl","4_SLK1", "5_SLK2","2_2SLT"), ] #Selects the rows with the specified sample IDs, effectively filtering the samples data frame to include only these particular samples for further analysis.

txi_subset4 <- txi
txi_subset4$counts <- txi$counts[, rownames(samples_subset4)]
txi_subset4$abundance <- txi$abundance[, rownames(samples_subset4)]
txi_subset4$length <- txi$length[, rownames(samples_subset4)]

dds_subset4 <- DESeqDataSetFromTximport(txi_subset4, colData = samples_subset4, design = ~ condition)
dds_subset4 <- DESeq(dds_subset4)
res_subset4 <- results(dds_subset4, contrast = c("condition", "Treatment", "Control"))
res_subset4 <- as.data.frame(res_subset4)
res_subset4$gene_id <- rownames(res_subset4)
res_with_metadata_subset4 <- merge(res_subset4, gene_metadata, by = "gene_id", all.x = TRUE)
write.csv(res_with_metadata_subset4, file = "~/Documents/DESeq2_results_subset4.csv", row.names = FALSE)




