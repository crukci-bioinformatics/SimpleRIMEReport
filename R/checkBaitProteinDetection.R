#' Set the possible negative control sample names
getNegCtrls <- function() {
    c("igg", "control", "empty", "none")
}

#' Check that the bait protein is present in the samples protein table
#'
#' This function checks that the bait protein Uniprot ID is present in the
#' protein table's `Accession` column
#' @name checkProtein
#' @import stringr
checkProtein <- function(SampleName, ProteinFile, BaitProteinUniProtID, ...) {
    proteinAcc <- read_tsv(ProteinFile, col_type = cols()) %>%
        pull(Accession)

    if (!BaitProteinUniProtID %in% proteinAcc) {
        message("The bait protein ",
                BaitProteinUniProtID,
                " was not detected in the protein table for sample ",
                SampleName,
                ".")
    }
}

#' Check that the bait protein is present in the samples peptide table
#'
#' This function checks that the bait protein Uniprot ID is present in the
#' peptide table's `Master Protein Accessions` column
#' @name checkPeptide
#' @import stringr
checkPeptide <- function(SampleName, PeptideFile, BaitProteinUniProtID, ...) {
    peptideAcc <- read_tsv(PeptideFile, col_type = cols()) %>%
        pull(`Master Protein Accessions`)

    if (!BaitProteinUniProtID %in% peptideAcc) {
        message("The bait protein ",
                BaitProteinUniProtID,
                " was not detected in the peptide table for sample ",
                SampleName,
                ".")
    }
}

#' Filter out rows for samples where the bait protein was not detected
#'
#' This function loads all the protein data for each sample and then checks that
#' for all samples the bait protein is detected in both tables.
#' @name checkBaitDetection
#' @import purrr
#' @import dplyr
#' @importFrom magrittr %>%
checkBaitDetection <- function(sampleTab) {
    negCtrls <- getNegCtrls()
    sampleTab %>%
        filter(!BaitProteinName %in% negCtrls) %>%
        pwalk(checkProtein) %>%
        pwalk(checkPeptide)
}
