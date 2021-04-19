
getProteinTables <- function(sampleTab, dataDir, annot){
    loadProteinData(sampleTab, dataDir) %>%   
        map(left_join, annot, by=c("Accession"="Accessions")) %>% 
        map(select, Accession, Gene, `Gene Symbol`=GeneSymbol, 
            Description, everything())
}