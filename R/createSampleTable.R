#' Extract information from file name
#'
#'  This function splits the file name and extracts the sample name and protein
#'
#' @import dplyr
#' @import stringr
#' @import tibble
#' @importFrom magrittr %>%
filenametotable <- function(fileName) {
    dat <- str_split(fileName, "_")[[1]]
    len <- length(dat)
    tibble(ProteinFile = fileName,
           PeptideFile = str_replace_all(ProteinFile,
                                         c(Proteins.txt = "PeptideGroups.txt",
                                           proteingroups.txt = "psms.txt")),
           SampleName = str_c(dat[6:(len-2)], collapse = "_"),
           BaitProteinName = dat[len - 1])
}

#' Create the sample table based on files in the data directory
#'
#' This function scans the contents of a directory for the presence of Proteome
#' Discoverer results files. It will detect both PD version 1.4 and PD version
#' 2.4 files.
#'
#' @param dataDir character; Path to directory containing the Proteome
#'   Discoverer output files
#' @details The files required for each sample are a Protein table and a Peptide
#' table. For PD version 1.4 filenames are expected to have the following naming
#' convention: '
#' 
#' QE_HF_<Date>_<PR_number>_<Initials>_<SampleName>_<BaitProteinName>_proteingroups.txt
#' QE_HF_<Date>_<PR_number>_<Initials>_<SampleName>_<BaitProteinName>_psms.txt
#'
#' For PD version 2.4 filenames are expected to have the following naming
#' convention:
#'
#' QE_HF_<Date>_<PR_number>_<Initials>_<SampleName>_<BaitProteinName>_Proteins.txt
#' QE_HF_<Date>_<PR_number>_<Initials>_<SampleName>_<BaitProteinName>_PeptideGroups.txt
#'
#' The function creates a new file, in the sample directory, that contains a
#' sample table, which is used to generate the report.
#'
#' The sample table will contain 5 columns:
##' \itemize{
##'  \item{ ProteinFile - the file name for 'protein groups' data }
##'  \item{ PeptideFile - the file name for 'psms' data }
##'  \item{ SampleName - the name for the sample - this is used for the plot
##'                      titles and the names of the worksheets }
##'  \item{ BaitProteinName - the name for the bait protein, e.g. "FOXA1",
##'                           or one of "IgG", "Control", or "None". }
##'  \item{ BaitProteinUniProtID - the UniProt ID for the bait protein, e.g.
##'                                P03372, P55317 this is used to look up the
##'                                fasta sequence }
##'  \item{ SampleType - "Control" or "Other" - this is used to determine which
##'                      samples are used to derive the "non-specific" proteins
##'                      for filtering.}
##' }
#'
#' The script will attempt to deduce the SampleName and BaitProteinName from the
#' the file name. It will also try to determine the UniProt ID from the bait
#' protein name. If the BaitProteinName is associated with multiple UniProt IDs,
#' the script will return all possibilities.
#'
#' It is important to check the details contained in the table, and modify it if
#' necessary, prior to generating the report. In particular check the bait
#' protein name and UniProt ID. The negative control should be called "Control",
#' "Empty" or "IgG" and have no UniProtID.
#' @return A sample table suitable for input into the `createReport` function. A
#'   copy of the table it also written to the data directory as a comma
#'   separated table.
#'
#' @examples
#'
#' @import dplyr
#' @import purrr
#' @import stringr
#' @import readr
#' @importFrom magrittr %>%
#'
#' @export createSampleTable
createSampleTable <- function(dataDir) {
    if (!dir.exists(dataDir)) {
        stop("The data directory provided - ", dataDir, " - cannot be found")
    }
    protFiles <- list.files(dataDir,
                            pattern = "_Proteins.txt$|_proteingroups.txt$")
    if (length(protFiles) == 0) {
        stop("No '_Protein.txt' or '_proteingroups.txt' files were found in ",
             dataDir)
    }
    tab <- map_df(protFiles, filenametotable)

    #Note `annot` is loaded from data/uniprot_annotation.rda
    data(uniprot_annotation)

    uid <- filter(annot, GeneSymbol %in% tab$BaitProteinName) %>%
        select(BaitProteinName = GeneSymbol, Accession) %>%
        group_by(BaitProteinName) %>%
        summarise(BaitProteinUniProtID = str_c(Accession, collapse = "/")) %>%
        ungroup()

    negCtrls <- getNegCtrls()
    tab <- left_join(tab, uid, by = "BaitProteinName") %>%
        mutate(SampleType = ifelse(str_to_lower(BaitProteinName) %in% negCtrls,
                                   "Control",
                                   "Other"))

    outnam <- str_c(dataDir, "/", basename(dataDir), "_sample_details.csv")
    write_csv(tab, outnam)
    message("A copy of the sample table has been saved at ", outnam)
    tab
}
