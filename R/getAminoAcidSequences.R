#' Create fasta table
#' 
#' Creates a table with the fasta sequence for each bait protein for the report
#' @name createFastaTable
#' 
#' @import dplyr
#' @import purrr
#' @import readr
#' @importFrom Biostrings AAStringSet
createFastaTable <- function(sampleTab){
    sampleTab %>% 
        select(BaitProteinUniProtID) %>% 
        distinct() %>% 
        filter(!is.na(BaitProteinUniProtID)) %>% 
        left_join(annot, by=c("BaitProteinUniProtID"="Accessions")) %>% 
        mutate(ProteinSeq = map(Sequence, AAStringSet)) %>% 
        select(BaitProteinUniProtID, ProteinSeq)
}