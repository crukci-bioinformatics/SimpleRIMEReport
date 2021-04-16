# split the file name and extract the sample name and protein
filenametotable <- function(fileName){
    dat <- str_split(fileName, "_")[[1]]
    len <- length(dat)
    tibble(ProteinFile=fileName,
           PeptideFile=str_replace(ProteinFile, "proteingroups", "psms"), 
           SampleName=str_c(dat[6:(len-2)], collapse="_"),
           BaitProteinName=dat[len-1])
}

#' Create the sample table based on files in the data directory
#'
#' This function scans the contents of a directory for the presence of Proteome
#' Discoverer results files. The files required for each sample are a Protein
#' table and a Peptide table. The filenames are expected to have the following
#' naming convention:
#'
#' QE_HF_<Date>_<PR_number>_<Initials>_<SampleName>_<BaitProteinName>_Proteins.txt
#' QE_HF_<Date>_<PR_number>_<Initials>_<SampleName>_<BaitProteinName>_PeptideGroups.txt
#'
#' The function creates a new file, in the sample directory, that contains a
#' sample table, which is used to generate the report. 
#' 
#' The sample table will contain 5 columns:
#' ProteinFile - the file name for 'protiengroups' data
#' PeptideFile - the file name for 'psms' data
#' SampleName - the name for the sample - this is used for the plot titles and
#'              the names of the worksheets
#' BaitProteinName - the protein name for the bait protein, e.g. ER, FOXA1 etc. 
#'                   This is only used in the plot titles, it could be anything
#' BaitProteinUniProtID - the UniProt ID for the bait protein, e.g. P03372, P55317
#'                       this is used to look up the fasta sequence
#' The script will attempt to guess the SampleName and BaitProteinName from the
#' the protein file name.
#' 
#' It is important to check the details contained in the table, and modify if
#' necessary, prior to generating the report. In particular check the bait
#' protein name and UniProt ID. The negative control should be called "Control"
#' or "IgG" and have no UniProtID.
#'
#' @param dataDir character; Path to directory containing the Proteome
#'   Discoverer output files
#' @return A sample table suitable for input into the `createReport` function. A
#'   copy of the table it also written to the data directory as a comma 
#'   separated table.
#' @examples
#'
#'
#' @import tidyverse
#'
#' @export createSampleTable
createSampleTable <- function(dataDir){
    tab <- list.files(dataDir, pattern="proteingroups.txt$") %>%
        map_df(filenametotable) %>% 
        mutate(BaitProteinUniProtID="")
    outnam <- str_c(dataDir, "/", basename(dataDir), "_sample_details.csv")
    write_csv(tab, outnam)
    return(tab)
}
