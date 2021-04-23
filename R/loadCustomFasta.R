#' Add data from custom fasta to annotation
#' 
#' This function reads annotations and sequences from fasta files and adds them
#' to the annotation used for the report.
#' The Protein accession, Gene ID and description are taken from the name for 
#' each sequence. the expected format is
#'    >sp|<Accession>|<GeneID> Description
#' e.g.
#'    >sp|NUP98-HOXA9|NUP98-HOXA9_HUMAN OS=Homo sapiens
#' where the Accession is NUP98-HOXA9 and the Gene ID is NUP98-HOXA9_HUMAN. Any
#' further information is add to the Description.
#' @name loadCustomFasta
#' @importFrom Biostrings readAAStringSet
#' @import dplyr
#' @import tibble
#' @import tidyr
loadCustomFasta <- function(additionalFasta){
    fasta <- readAAStringSet(additionalFasta)
    
    addAnnot <- tibble(N = names(fasta), Sequence = as.data.frame(fasta)$x) %>% 
        separate(N, 
                 into=c("DB", "Accession", "Gene", "Description"), 
                 sep=" |\\|",
                 extra = "merge") %>% 
        select(-DB)

    annotTab <- bind_rows(annot, addAnnot)
    return(annotTab)
} 
