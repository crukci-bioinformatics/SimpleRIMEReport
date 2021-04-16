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

getProteinTables <- function(sampleTab, dataDir, annot){
    loadProteinData(sampleTab, dataDir) %>%   
        map(left_join, annot, by=c("Accession"="Accessions")) %>% 
        map(select, Accession, Gene, `Gene Symbol`=GeneSymbol, 
            Description, everything())
}

# load peptide data. Each table into a element of a list
getPeptideTables <- function(sampleTab, dataDir){
    str_c(dataDir, "/", sampleTab$PeptideFile) %>% 
        map(read_tsv, col_types = cols()) %>% 
        set_names(sampleTab$SampleName)
}

