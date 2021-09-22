#' Check that the control samples do not have a UniProt ID associated with them
#' 
#' This function checks that the bait protein Uniprot ID column is empty (NA) 
#' for control samples
#' @name check
#' @import stringr
#' @import dplyr
checkControlIDs <- function(sampleTable){
    negCtrls <- c("igg", "control", "empty")
    allNA <- sampleTable %>%
        filter(str_to_lower(BaitProteinName)%in%negCtrls) %>%
        pull(BaitProteinUniProtID) %>% 
        is.na() %>% 
        all()
    if(!allNA){
        stop("There are control samples (IgG/Control/Emtpy) that have bait ",
             "protein Uniprot IDs. Please check your sample table.")
    }
}
