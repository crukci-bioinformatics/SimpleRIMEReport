#' Check that the bait protein is present in the samples protein table
#' 
#' This function checks that the bait protein Uniprot ID is present in the 
#' protein and peptide tables `Accession` column, or that the sample is a 
#' negative control.
#' @name checkProtein
#' @import stringr
checkProtein <- function(SampleName, BaitProteinName, BaitProteinUniProtID, 
                          ProteinDat, PeptideDat, ...){
    negCtrls <- c("igg", "control", "empty")
    isGoodprot <- str_to_lower(BaitProteinName)%in%negCtrls |
                BaitProteinUniProtID%in%ProteinDat$Accession
    if(!isGoodprot){
        message("The bait protein ", BaitProteinName, ":", BaitProteinUniProtID,
                " was not detected in the protein table for sample ", SampleName,
                ".") 
    }
    isGoodpep <- str_to_lower(BaitProteinName)%in%negCtrls |
        BaitProteinUniProtID%in%PeptideDat$`Master Protein Accessions`
    if(!isGoodpep){
        message("The bait protein ", BaitProteinName, ":", BaitProteinUniProtID,
                " was not detected in the peptide table for sample ", SampleName,
                ".") 
    }
    return(isGoodprot & isGoodpep)
}

#' Filter out rows for samples where the bait protein was not detected
#' 
#' This function loads all the protein data for each sample and then checks that
#' for all samples the bait protein is detected in both tables.
#' @name checkBaitDetection
#' @import purrr
#' @import dplyr
#' @importFrom magrittr %>%
checkBaitDetection <- function(sampleTab){
    areGood <- sampleTab %>% 
        mutate(ProteinDat = map(ProteinFile, read_tsv, col_type = cols())) %>%
        mutate(PeptideDat = map(PeptideFile, read_tsv, col_type = cols())) %>%
        pmap(checkProtein) %>%   
        unlist()
    if(!all(areGood)){
        stop("There are samples for which the bait protein provided was not ",
             "detected in the data. Please check your sample table.")
    }
}
