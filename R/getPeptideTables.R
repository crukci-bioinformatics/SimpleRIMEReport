#' Get peptide tables
#' 
#' This function reads each Peptide data file and creates list of peptide level 
#' data tables -  one for each bait and the control.
#' @name getPeptideTables
#' @import readr
#' @import purrr
getPeptideTables <- function(sampleTab){
    map(sampleTab$PeptideFile, read_tsv, col_types = cols()) %>% 
        set_names(sampleTab$SampleName)
}