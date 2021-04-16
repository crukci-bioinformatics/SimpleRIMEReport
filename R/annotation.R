loadAnnotation <- function(sampleTab, dataDir, species_id, upCon){
    columns <- c("ENTRY-NAME", "PROTEIN-NAMES", "GENES")
    loadProteinData(sampleTab, dataDir) %>%   
        map(pull, var=1) %>%
        unlist() %>%
        unique() %>%
        UniProt.ws::select(upCon, ., columns, "UNIPROTKB") %>%
        as_tibble() %>% 
        mutate(GeneSymbol = gsub(" .*", "", GENES)) %>% 
        select(Accessions = "UNIPROTKB", Gene = "ENTRY-NAME",
               Description = "PROTEIN-NAMES", GeneSymbol)
}

makeAnnotation <- function(sampleTab, dataDir){
    str_c(dataDir, "/", sampleTab$ProteinFile) %>% 
        map_df(read_tsv, col_types = cols()) %>%  
        select(Accessions=Accession) %>%  
        distinct() %>%   
        left_join(readRDS("annotation/SwissProtAnnotationTable.rds"), 
                  by="Accessions") %>% 
        select(-Sequence)
}
