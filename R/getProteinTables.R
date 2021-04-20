#' Get protein tables
#' 
#' This function reads each Protein data file and creates list of protein level 
#' data tables with annotations -  one for each bait and the control.
#' @name getProteinTables
#' @import dplyr
#' @import readr
#' @import purrr
getProteinTables <- function(sampleTab){
    loadProteinData(sampleTab) %>%   
        map(left_join, annot, by=c("Accession"="Accessions")) %>% 
        map(select, Accession, Gene, `Gene Symbol`=GeneSymbol, 
            Description, everything(), -Sequence)
}