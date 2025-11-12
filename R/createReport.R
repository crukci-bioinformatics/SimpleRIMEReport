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
                         outputFileName = NULL) {

    sampleTable <- sampleTable %>%
        mutate(across(ends_with("File"), ~str_c(dataDir, "/", .x)))

    #Note `annot` is loaded from data/uniprot_annotation.rda
    data(uniprot_annotation)
    annotTab <- annot

    # Load  any additional fasta files and add to the annotation
    if (!is.null(additionalFasta)) {
        annotTab <- addAnnotation(annotTab, additionalFasta)
    }

    # Check that the annotation has the baits
    checkAnnForBait(sampleTable, annotTab)

    # load peptide data and check for bait proteins
    sampleTable <- getPeptideTables(sampleTable)
    checkBaitDetection(sampleTable, type = "peptide")

    # load the protein data and check for bait proteins
    sampleTable <- getProteinTables(sampleTable, annotTab)
    checkBaitDetection(sampleTable, type = "protein")

    # filter the protein data and check for bait proteins
    sampleTable <- filterProteinTables(sampleTable)
    checkBaitDetection(sampleTable, type = "filtered protein")

    # make merged tables
    combined_data <- combineProteinTables(sampleTable)
    filtered_data <- filterCombinedTable(sampleTable)

    # make plots
    coveragePlots <- makeCoveragePlots(sampleTable, annotTab)

    # make protein data list
    proteinData <- sampleTable$ProteinTable
    names(proteinData) <- sampleTable$SampleName

    # make workbook
    if (is.null(outputFileName)) {
        outputFileName <- str_c(dataDir,
                                "/",
                                basename(dataDir),
                                ".Simple_RIME_report.xlsx")
    }
    makeWorkBook(outputFileName,
                 plotsList = coveragePlots,
                 protData = proteinData,
                 combData = combined_data,
                 filtData = filtered_data)
    message("Report written to ", outputFileName)
}