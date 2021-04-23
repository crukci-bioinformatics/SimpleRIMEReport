#' Create the Simple RIME report
#'
#' This function generates report as an Excel workbook. The input is a sample
#' table.
#'
#' @param sampleTable data.frame; the sample table
#' @param dataDir character; Path to directory containing the Proteome
#'   Discoverer output files
#' @param additionalFasta character; a vector of paths to custom fasta files to
#'   be included in the annotation (optional)
#' @param outputFileName charcter; an output file name (optional)
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

createReport <- function(sampleTable, 
                         dataDir, 
                         additionalFasta = NULL,
                         outputFileName = NULL){
    
    sampleTable <- sampleTable %>% 
        mutate(across(ends_with("File"), ~str_c(dataDir, "/", .x)))

    # check that all the bait proteins were detected
    sampleTable <- checkBaitDetection(sampleTable)

    if(all(is.na(sampleTable$BaitProteinUniProtID))){
        stop("There are no samples in which the bait protein provided was ",
             "detected.")
    }

    # Load  any additional fasta files and add to the annotation
    if(!is.null(additionalFasta)){ 
        annotTab <- loadCustomFasta(additionalFasta) 
    }else{
        annotTab <- annot
        }
    # load data
    protein_data <- getProteinTables(sampleTable, annotTab)
    peptide_data <- getPeptideTables(sampleTable)

    # make merged tables
    combined_data <- combineProteinTables(sampleTable, annotTab)
    filtered_data <- filterCombinedTable(sampleTable, annotTab)

    # make plots
    coveragePlots <- makeCoveragePlots(sampleTable, peptide_data, annotTab)

    # make workbook
    if(is.null(outputFileName)){
        outputFileName <- str_c(dataDir, 
                                "/", 
                                basename(dataDir),
                                ".Simple_RIME_report.xlsx")
    }
    makeWorkBook(outputFileName,
                 coveragePlots,
                 protein_data,
                 combined_data,
                 filtered_data)
    message("Report written to ", outputFileName)
}
