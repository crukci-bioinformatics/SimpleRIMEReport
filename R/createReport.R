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
#' @import GenomicRanges
#' @import Biostrings
#' @import openxlsx
#' @import cowplot
#' @import tidyverse
#'
#' @export createReport

createReport <- function(sampleTable, dataDir){

    # Load from preconstructed SwissProt file
    annotTab <- makeAnnotation(sampleTable, dataDir)
    
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