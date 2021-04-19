#' Create annotation table
#' 
#' Creates an annotation table for the report
#' @name makeAnnotation
#' 
#' @import dplyr
#' @import purrr
#' @import readr

makeAnnotation <- function(sampleTab, dataDir){
    map_df(sampleTab$ProteinFile, read_tsv, col_types = cols()) %>%
        select(Accessions=Accession) %>%  
        distinct() %>%   
        left_join(annot, by="Accessions") %>%
        select(-Sequence)
}
