#' Add sample name to column header
#'
#' Add the sample name to all column headers except "Accession"
#' @name getNonSpecificProteins
#' @import dplyr
#' @import stringr
addSampleName <- function(tab, sampleName) {
    rename_with(tab,
                ~str_c(sampleName, " ", .),
                c("Coverage", "Unique Peptides"))
}

#' Combine all protein data into a single table
#'
#' Read all the protein data tables, combine them into a single table - add the
#' sample name to the beginning of each header.
#' @name combineProteinTables
#' @import dplyr
#' @import purrr
combineProteinTables <- function(sampleTab) {
    pull(sampleTab, ProteinTable) %>%
        map2(sampleTab$SampleName, addSampleName) %>%
        reduce(full_join, by = c("Accession",
                                 "Gene",
                                 "Gene Symbol",
                                 "Description"))
}

#' Get Non-specific proteins
#'
#' Read all protein tables for negative controls and compile a vector of
#' the proteins detected. These will be considered as present due to
#' non-specific binding and removed from the final filtered protein table.
#' @name getNonSpecificProteins
#' @import dplyr
#' @import stringr
#' @import purrr
getNonSpecificProteins <- function(sampleTab) {
    sampleTab %>%
        filter(SampleType == "Control") %>%
        pull(ProteinTable) %>%
        reduce(full_join, by = "Accession") %>%
        pull(Accession) %>%
        unique()
}

#' Create filtered protein table
#'
#' Read all the protein data tables for the baited samples, combine them into a
#' single table. Then filter out all proteins that are detected in the negative
#' controls.
#' @name filterCombinedTable
#' @import dplyr
#' @import stringr
#' @import purrr
filterCombinedTable <- function(sampleTab) {
    non_specific <- getNonSpecificProteins(sampleTab)

    sampleTab %>%
        filter(SampleType != "Control") %>%
        combineProteinTables() %>%
        filter(!Accession %in% non_specific)
}