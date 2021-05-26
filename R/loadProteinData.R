#' Select columns from Protein data table
#' 
#' This function selects the columns to keep from the PD protein table.
#' @name selectColumns
#' @import dplyr
#' @import tidyselect
selectColumns <- function(x){
    dplyr::select(x, 
                  Accession,
                  Coverage = contains("Coverage"),
                  `Unique Peptides`=`# Unique Peptides`,
                  any_of("Master"))
}

#' Filter Protein data table for MasterProteins
#'
#' This function filters on the  columns "Master" in PD v2.4 files to keep only
#' those rows that are marked "IsMasterProtein". It then removes the "Master"
#' column as it is no longer needed.
#' @name filterMasterProteins
#' @import dplyr
#' @import tidyselect
filterMasterProteins <- function(x){
    dplyr::filter(x, across(any_of("Master"), ~.x == "IsMasterProtein")) %>%
        dplyr::select(-any_of("Master"))
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
loadProteinData <- function(sampleTab){
    sampleTab$ProteinFile %>% 
        map(read_tsv, col_types = cols()) %>% 
        map(selectColumns) %>% 
        map(filterMasterProteins) %>%
        set_names(sampleTab$SampleName)
}

