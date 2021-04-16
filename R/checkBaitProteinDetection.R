#' Check that the bait protein is present in the samples protein table
#' 
#' This function checks that the bait protein Uniprot ID is present in the 
#' protein tables `Accession` column, or that the sample is a negative control.
#' @name checkProtein
#' @import stringr
checkProtein <- function(SampleName, BaitProteinName, BaitProteinUniProtID, 
                          ProteinDat, ...){
    negCtrls <- c("igg", "control", "empty")
    isGood <- str_to_lower(BaitProteinName)%in%negCtrls |
                BaitProteinUniProtID%in%ProteinDat$Accession
    if(!isGood){
        warning("The bait protein ", BaitProteinName, " was not detected in ",
                "the sample ", SampleName, ". The results for this sample ",
                "will be excluded from the report.") 
    }
    return(isGood)
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
        mutate(ProteinDat = map(ProteinFile, read_tsv)) %>%   
        #dplyr::select(SampleName, BaitProteinName, BaitProteinUniProtID, ProteinDat) %>% 
        pmap(checkProtein) %>%   
        unlist()
    return(sampleTab[areGood,])
}
