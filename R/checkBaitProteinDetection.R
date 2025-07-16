#' Set the possible negative control sample names
getNegCtrls <- function() {
    c("igg", "control", "empty", "none")
}

#' Check that the bait protein is present in the table
#'
#' This function checks that the bait protein Uniprot ID is present in the
#' table's `Accessions` column
#' @name checkBait
#' @import stringr
checkBait <- function(SampleName, Table, BaitProteinUniProtID, Type, ...) {
    protAcc <- Table$Accession

    if (!BaitProteinUniProtID %in% protAcc) {
        message("The bait protein ",
                BaitProteinUniProtID,
                " was not detected in the ",
                Type,
                " table for sample ",
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
checkBaitDetection <- function(sampleTab, type = "protein") {
    negCtrls <- getNegCtrls()
    tabCol <- ifelse(type == "peptide", "PeptideTable", "ProteinTable")
    sampleTab %>%
        mutate(across(BaitProteinName, str_to_lower)) %>%
        filter(!BaitProteinName %in% negCtrls) %>%
        rename(Table = !!tabCol) %>%
        mutate(Type = type) %>%
        pwalk(checkBait)
}
