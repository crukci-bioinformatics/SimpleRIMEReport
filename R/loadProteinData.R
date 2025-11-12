#' Select columns from Protein data table
#'
#' This function selects the columns to keep from the PD protein table.
#' @name selectColumns
#' @import dplyr
#' @import tidyselect
selectColumns <- function(x) {
    dplyr::select(x,
                  Accession,
                  `Unique Peptides` = `# Protein Unique Peptides`,
                  any_of("Master"))
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
loadProteinTable <- function(filePath, annotTab) {
    read_tsv(filePath, show_col_types = FALSE) %>%
        selectColumns() %>%
        left_join(annotTab, by = c("Accession" = "Accession")) %>%
        select(Accession,
               Gene,
               `Gene Symbol` = GeneSymbol,
               Description,
               everything(),
               -Sequence)
}

#' Get protein tables
#'
#' This function reads each Protein data file and creates list of protein level
#' data tables with annotations -  one for each bait and the control.
#' @name getProteinTables
#' @import dplyr
#' @import readr
#' @import purrr
getProteinTables <- function(sampleTab, annotTab) {
    sampleTab %>%
        mutate(ProteinTable = map(ProteinFile,
                                  ~loadProteinTable(.x, annotTab = annotTab)))
}

#' Filter Protein data table for MasterProteins
#'
#' This function filters on the  columns "Master" in PD v2.4 files to keep only
#' those rows that are marked "IsMasterProtein". It then removes the "Master"
#' column as it is no longer needed.
#' @name filterMasterProteins
#' @import dplyr
#' @import tidyselect
filterProteins <- function(x) {
    dplyr::filter(x, if_any(contains("Master"), ~.x == "IsMasterProtein")) %>%
        dplyr::select(-any_of("Master"))
}

#' Filter the protein tables
#'#'
#' This function filters the protein tables in the sample table
#' @name filterProteinTables
#' @import purrr
filterProteinTables <- function(sampleTab) {
    sampleTab %>%
        mutate(ProteinTable = map(ProteinTable, filterProteins))
}
