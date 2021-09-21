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
        message("The bait protein ", BaitProteinName, " was not detected in ",
                "the protein table for sample ", SampleName, ". The results ",
                "for this sample will be excluded from the report.") 
    }
    isGoodpep <- str_to_lower(BaitProteinName)%in%negCtrls |
        BaitProteinUniProtID%in%PeptideDat$`Master Protein Accessions`
    if(!isGoodpep){
        message("The bait protein ", BaitProteinName, " was not detected in ",
                "the peptide table for sample ", SampleName, ". The results ",
                "for this sample will be excluded from the report.") 
    }
    return(isGoodprot & isGoodpep)
}

#' Filter out rows for samples where the bait protein was not detected
#' 
#' This function loads all the protein data for each sample and then eliminates
#' rows where the bait protein was not detected.
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
    return(sampleTab[areGood,])
}
