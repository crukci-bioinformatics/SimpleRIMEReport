#' Add data from custom fasta to annotation
#'
#' This function will add the annotation from custom fasta file to the
#' annotation used for the report if necessary. The Protein accession, Gene ID
#' and description are taken from the name for each sequence. The default
#' expected format is:
#'    >{Database}@{Accession}@{GeneID}@{Description}
#' Where the separators (here @) can be "|" or " ".
#' e.g.
#'    >sp|NUP98-HOXA9|NUP98-HOXA9_HUMAN OS = Homo sapiens
#' where the Accession is NUP98-HOXA9 and the Gene ID is NUP98-HOXA9_HUMAN. Any
#' further information is add to the Description.
#' If an accession in the custom fasta is also present in the original
#' annotation, it will be removed from the original annotation on the basis that
#' the (presumably) modified entry is the one that is required.
#' This is the default header for protein sequence fasta files from Uniprot
#' @name addAnnotation
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import tidyr
addAnnotation <- function(annotTab, additionalFasta) {
    fasta <- readAAStringSet(additionalFasta)

    dbCols <- c("DB", "Accession", "Gene", "Description")
    addAnnot <- tibble(N = names(fasta), Sequence = as.data.frame(fasta)$x) %>%
        separate_wider_delim(N,
                             delim = regex(" |\\|"),
                             names = dbCols,
                             too_few = "align_start",
                             too_many = "merge") %>%
        select(-DB)

    # Check the the additional annotation has an Accessions column
    msg <- vector()
    if (!"Accession" %in% colnames(addAnnot)) {
        msg <- str_c("The additional annotation does not have an 'Accession'",
                     " column. \n")
    }
    # Check that there are no missing accessions
    if (any(is.na(addAnnot$Accession))) {
        msg <- c(msg,
                 "The additional annotation has missing accessions. \n")
    }
    # Check that the additional annotation has unique accessions
    if (any(duplicated(addAnnot$Accession))) {
        msg <- c(msg,
                 "The additional annotation contains duplicated accessions. \n")
    }
    if (length(msg) > 0) {
        write_tsv(addAnnot, "Annotation_with_errors.tsv")
        stop(msg,
             "Please check the format of the fasta file headers ",
             "(see vignette for details).\n",
             "The derived annotation has been saved to ",
             "'Annotation_with_errors.tsv'.")
    }

    # remove duplicated accessions from the original annotation
    annotTab <- annotTab %>%
        filter(!Accession %in% addAnnot$Accession)

    bind_rows(annotTab, addAnnot)
}

#' Check bait in annotation
#' This function checks that all the baits in the sample table are present in
#' the annotation. If there are any missing baits the function throws and error
#' and exits the script.
#' @name checkAnnForBait
checkAnnForBait <- function(sampleTab, annotTab) {
    baitIDs <- sampleTab %>%
        filter(!is.na(BaitProteinUniProtID)) %>%
        pull(BaitProteinUniProtID)
    if (!all(baitIDs %in% annotTab$Accession)) {
        missingBaits <- baitIDs[!(baitIDs %in% annot$Accession)] %>%
            unique() %>%
            str_c("\t", .) %>%
            str_c(collapse = "\n")
        warning("The following the bait IDs are not found in the annotation:\n",
                missingBaits,
                "\nPerhaps you need to provide a custom fasta for annotation.")
    }
}
