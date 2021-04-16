# get fasta sequence from uniprot using Accession
getFasta <- function(UniProtAccession, species_id, upCon){
    #hs <- UniProt.ws::UniProt.ws(taxId = species_id)
    aaSeq <- UniProt.ws::select(upCon, UniProtAccession, "SEQUENCE", "UNIPROTKB") %>%
        pull(SEQUENCE) %>%   
        AAStringSet()
}

makeFastaTable <- function(sampleTab, species_id, upCon){
    sampleTab %>% 
        select(BaitProteinUniProtID) %>% 
        distinct() %>% 
        filter(!is.na(BaitProteinUniProtID)) %>% 
        mutate(ProteinSeq=map(BaitProteinUniProtID, 
                              getFasta, 
                              species_id=species_id,
                              upCon=upCon))
}

createFastaTable <- function(sampleTab){
    sampleTab %>% 
        select(BaitProteinUniProtID) %>% 
        distinct() %>% 
        filter(!is.na(BaitProteinUniProtID)) %>% 
        filter(!str_detect(BaitProteinUniProtID, "^[Ii][Gg][Gg]$")) %>% 
        left_join(readRDS("annotation/SwissProtAnnotationTable.rds"), 
                  by=c("BaitProteinUniProtID"="Accessions")) %>% 
        mutate(ProteinSeq=map(Sequence, AAStringSet)) %>%  
        select(BaitProteinUniProtID, ProteinSeq)
}