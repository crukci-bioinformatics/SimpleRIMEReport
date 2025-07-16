#' Load pepti
loadPeptideTable <- function(filePath) {
    # Read the peptide data file
    read_tsv(filePath, show_col_types = FALSE) %>%
        rename_with(~str_replace(.x,
                                 pattern = ".*Accessions$",
                                 replacement = "Accession"))
}

#' Get peptide tables
#'
#' This function reads each Peptide data file into a list column in the sample
#table. 
#' @name getPeptideTables
#' @import readr
#' @import purrr
getPeptideTables <- function(sampleTab){
    sampleTab %>%
        mutate(PeptideTable = map(PeptideFile, loadPeptideTable))
}
