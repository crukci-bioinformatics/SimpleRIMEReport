# add sample name to column headers
addSampleName <- function(tab, sampleName){
    tab %>% 
        rename_at(vars(-Accession), 
                  list(~str_c(sampleName, " ", .)))
}

combineProteinTables <- function(sampleTab, dataDir, annot){
    loadProteinData(sampleTab, dataDir) %>%   
        map2(sampleTable$SampleName, addSampleName) %>% 
        reduce(full_join, by="Accession") %>% 
        left_join(annot, by=c("Accession"="Accessions")) %>% 
        select(Accession, Gene, `Gene Symbol`=GeneSymbol, Description, everything())
}

filterCombinedTable <- function(sampleTab, dataDir, annot){
    non_specific <- sampleTab %>% 
        filter(str_detect(BaitProteinName, "^[Ii][Gg][Gg]$")) %>%  
        loadProteinData(dataDir) %>%   
        reduce(full_join, by="Accession") %>% 
        pull(Accession) %>% 
        unique()
    
    sampleTab <- sampleTab %>% 
        filter(!str_detect(BaitProteinName, "^[Ii][Gg][Gg]$"))
    loadProteinData(sampleTab, dataDir) %>%   
        map2(sampleTab$SampleName, addSampleName) %>% 
        reduce(full_join, by="Accession") %>% 
        left_join(annot, by=c("Accession"="Accessions")) %>% 
        select(Accession, Gene, `Gene Symbol`=GeneSymbol, Description, everything()) %>% 
        filter(!Accession%in%non_specific)
}