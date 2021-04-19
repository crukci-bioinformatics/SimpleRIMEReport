#' Select columns from Protein data table
#' 
#' This function selects the columns to keep from the PD protein table.
#' @name selectColumns
#' @import dplyr
selectColumns <- function(x){
    dplyr::select(x, 
                  Accession,
                  Coverage = `Coverage [%]`,
                  `Unique Peptides`=`# Unique Peptides`)
}

#' Load protein data table
#' 
#' Load protein data. Each table is read into an element of a list, which is 
#' named according to the sample names.
#' @name loadProteinData
#' @import stringr
#' @import purrr
#' @import dplyr
#' @import readr
#' @importFrom magrittr %>%
loadProteinData <- function(sampleTab, dataDir){
    sampleTab$ProteinFile %>% 
        map(read_tsv, col_types = cols()) %>% 
        map(selectColumns) %>% 
        set_names(sampleTab$SampleName)
}

