#' Check that the bait protein is present in the samples protein table
#' 
#' This function checks that the bait protein is present in the protein table
#' @name checkProteins
#' @import tidyverse
checkProteins <- function(SampleName, BaitProteinName, BaitProteinUniProtID, 
                          ProteinDat){
    isGood <- str_detect(BaitProteinName, "^[Ii][Gg][Gg]$") |
        BaitProteinUniProtID%in%ProteinDat$Accession
    if(!isGood){
        warning("The bait protein ", BaitProteinName, " was not detected in ",
                "the sample ", SampleName, ". The results for this sample ",
                "will be excluded from the report.") 
    }
    isGood
}

#' Filter out rows for samples where the bait protein was not detected
#' 
#' This function loads all the protein data for each sample and then eliminates
#' rows where the bait protein was not detected.
#' @name checkBaitDetection
#' @import tidyverse
checkBaitDetection <- function(sampleTab, dataDir){
    ProtData <- loadProteinData(sampleTab, dataDir)
    areGood <- sampleTab %>% 
        dplyr::select(SampleName, BaitProteinName, BaitProteinUniProtID) %>% 
        mutate(ProteinDat=ProtData) %>%   
        pmap(checkProteins) %>%   
        unlist()
    sampleTab[areGood,]
}