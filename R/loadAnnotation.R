#' Add data from custom fasta to annotation
#'
#' This function will add the annotation from custom fasta file to the
#' annotation used for the report if necessary. The Protein accession, Gene ID
#' and description are taken from the name for each sequence. the expected
#' format is:
#'    >sp|<Accession>|<GeneID> <Description>
#' e.g.
#'    >sp|NUP98-HOXA9|NUP98-HOXA9_HUMAN OS = Homo sapiens
#' where the Accession is NUP98-HOXA9 and the Gene ID is NUP98-HOXA9_HUMAN. Any
#' further information is add to the Description.
#' The function then checks that all the bait protein IDs are present in the
#' annotation Accessions
#' @name loadAnnotation
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import tidyr
loadAnnotation <- function(sampleTab, additionalFasta){
    #Note `annot` is loaded automatically from R/sysdata.rda

    if(!is.null(additionalFasta)){
        fasta <- readAAStringSet(additionalFasta)

        addAnnot <- tibble(N = names(fasta), Sequence = as.data.frame(fasta)$x) %>%
            separate(N,
                     into = c("DB", "Accession", "Gene", "Description"),
                     sep = " |\\|",
                     extra = "merge") %>%
            select(-DB)

        annot <- bind_rows(annot, addAnnot)
    }

    baitIDs <- sampleTab %>%
        filter(!is.na(BaitProteinUniProtID)) %>%
        pull(BaitProteinUniProtID)
    if(!all(baitIDs %in% annot$Accession)){
        missingBaits <- baitIDs[!(baitIDs %in% annot$Accession)] %>%
            unique() %>%
            str_c(collapse = "\n")
        stop("The following the bait IDs are not found in the annotation:\n",
             missingBaits,
             "\nPerhaps you need to provide a custom fasta for annotation.")
    }

    return(annot)
}
