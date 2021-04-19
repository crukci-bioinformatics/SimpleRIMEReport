#' Create the Simple RIME report
#'
#' This function generates report as an Excel workbook. The input is a sample
#' table.
#'
#' @param sampleTable data.frame; the sample table
#' @param dataDir character; Path to directory containing the Proteome
#'   Discoverer output files
#' @return NULL
#' @examples
#'
#' @import stringr
#'
#' @export createReport
 
# @import GenomicRanges
# @import Biostrings
# @import openxlsx
# @import cowplot

createReport <- function(sampleTable, dataDir){
    
    sampleTable <- sampleTable %>% 
        mutate(across(ends_with("File"), ~str_c(dataDir, "/", .x)))

    # check that all the bait proteins were detected
    sampleTable <- checkBaitDetection(sampleTable)

    if(all(is.na(sampleTable$BaitProteinUniProtID))){
        stop("There are no samples in which the bait protein provided was ",
             "detected.")
    }

    # Get annotation
    annotTab <- makeAnnotation(sampleTable)

    # Get fasta seq from SwissProt file
    fastaSeqTable <- createFastaTable(sampleTable)

    # load data
    protein_data <- getProteinTables(sampleTable, dataDir, annotTab)
    peptide_data <- getPeptideTables(sampleTable, dataDir)

    # make merged tables
    combined_data <- combineProteinTables(sampleTable, dataDir, annotTab)
    filtered_data <- filterCombinedTable(sampleTable, dataDir, annotTab)

    # make plots
    coveragePlots <- makeCoveragePlots(sampleTable, peptide_data, fastaSeqTable)

    # make workbook
    makeWorkBook(outputFileName,
                       coveragePlots,
                       protein_data,
                       combined_data,
                       filtered_data)
}
