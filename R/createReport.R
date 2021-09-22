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
#' @details The sample table should contain 5 columns:  
##' \itemize{
##'  \item{ ProteinFile - the file name for 'protiengroups' data }
##'  \item{ PeptideFile - the file name for 'psms' data }
##'  \item{ SampleName - the name for the sample - this is used for the plot
##'                      titles and the names of the worksheets }
##'  \item{ BaitProteinName - the protein name for the bait protein, e.g. ER,
##'                           FOXA1 etc. This is only used in the plot titles,
##'                           it could be anything }
##'  \item{ BaitProteinUniProtID - the UniProt ID for the bait protein, e.g.
##'                                P03372, P55317 this is used to look up the
##'                                fasta sequence }
##' }
#' The proteins are filtered on the "Master" column, if present, to only keep
#' proteins marked "IsMasterProtein".
#' 
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
    checkBaitDetection(sampleTable)

    # Check that the controls have not Bait Protein ID
    checkControlIDs(sampleTable)
    
    # Load  any additional fasta files and add to the annotation
    annotTab <- loadAnnotation(sampleTable, additionalFasta)
    
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
