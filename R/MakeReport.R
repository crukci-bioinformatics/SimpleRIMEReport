library(GenomicRanges)
library(Biostrings)
library(tidyverse)
library(openxlsx)
library(cowplot)
# Get functions
source("createSampleTable.R")
source("loadDataFunctions.R")
source("plotFunctions.R")
source("createWorkBookFunctions.R")

################################################################################
# 1. Set Input Parameters
################################################################################

data_directory <- "../TEST_RIME_SCRIPT/PR1001/"
# species_id <- as.numeric(9606)
# output file name
outputFileName <- str_c(data_directory, "/", 
                        basename(data_directory), ".report.xlsx")

################################################################################
# 2. Generate sample table
# This portion of the scripts uses the file names to generate a sample table
# Run the script and then open the table it saves in Excel and make any 
# changes or additions as needed.
# The table will be saved in the 'data_directory' as "XXXXX_sample_details.csv"
# The sample table will contain 5 columns:
# ProteinFile - the file name for 'protiengroups' data
# PeptideFile - the file name for 'psms' data
# SampleName - the name for the sample - this is used for the plot titles and
#              the names of the worksheets
# BaitProteinName - the protein name for the bait protein, e.g. ER, FOXA1 etc. 
#                   This is only used in the plot titles, it could be anything
# BaitProteinUniProtID - the UniProt ID for the bait protein, e.g. P03372, P55317
#                       this is used to look up the fasta sequence
# The script will attempt to guess the SampleName and BaitProteinName from the
# the protein file name.
# 
# The assumed structure of the protein file name is:
# QE_HF_<Date>_<PR_number>_<Initials>_<SampleName>_<BaitProteinName>_proteingroups.txt
# SampleName can include "_"
# e.g. 
#   QE_HF_04092018_PR1025_SAN_Sample1_VEH_ARID2_proteingroups.txt
# gives
#   SampleName="Sample1_VEH"
#   BaitProteinName="ARID2"
#
# You'll need to check the sample sheet and modify as necessary before 
# proceeding
################################################################################

createSampleTable(data_directory)

################################################################################
# 3. Generate the report
# Once you have adjusted the sample table run the rest of the code to generate
# the report
################################################################################

# get sample table
sampleTable <- str_c(data_directory, "/", 
                     basename(data_directory), "_sample_details.csv") %>% 
    read_csv(col_types = cols()) %>% 
    checkBaitDetection(data_directory)

if(all(str_detect(sampleTable$BaitProteinName, "^[Ii][Gg][Gg]$"))){
    stop("There are no samples in which the bait protein provided was detected.")
}

# make an annotation table.
# For each acc get annotations from UniProt
# uniProtConnection <- UniProt.ws(taxId = species_id)
# annotTab <- loadAnnotation(sampleTable, 
#                         data_directory, 
#                         species_id, 
#                         uniProtConnection)
# Get fasta seq from Uniprot
# fastaSeqTable <- makeFastaTable(sampleTable, species_id, uniProtConnection)

# make an annotation table.
# Load from preconstructed SwissProt file
annotTab <- makeAnnotation(sampleTable, data_directory)

# Get fasta seq from SwissProt file 
fastaSeqTable <- createFastaTable(sampleTable)

# load data
protein_data <- getProteinTables(sampleTable, data_directory, annotTab)
peptide_data <- getPeptideTables(sampleTable, data_directory)

# make merged tables
combined_data <- combineProteinTables(sampleTable, data_directory, annotTab)
filtered_data <- filterCombinedTable(sampleTable, data_directory, annotTab)

# make plots
coveragePlots <- makeCoveragePlots(sampleTable, peptide_data, fastaSeqTable)

# make workbook
makeWorkBook(outputFileName,
                   coveragePlots,
                   protein_data,
                   combined_data,
                   filtered_data)
