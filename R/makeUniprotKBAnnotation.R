library(tidyverse)

upURL <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release"
spURL <- str_c(upURL, "/knowledgebase/complete/uniprot_sprot.dat.gz") 

spTab <- read_tsv(spURL, col_names=FALSE, trim_ws = FALSE) %>%  
#spTab <- read_tsv("annotation/uniprot_sprot.dat", col_names=FALSE, trim_ws = FALSE) %>%  
    separate(X1, c("Code", "Value"), extra = "merge") %>%  
    filter(Code%in%c("ID", "AC", "GN", "DE", "SQ", "")) %>%  
    filter(Value!="") %>%  
    mutate(Gene=ifelse(Code=="ID", Value, NA)) %>%  
    mutate(Gene=str_remove(Gene, " .*")) %>%  
    fill(Gene, .direction = "down") %>%  
    filter(! (Code=="DE" & str_detect(Value, "RecName", negate = TRUE)))  %>%  
    mutate(Code=str_replace(Code, "^$", "SQ") ) %>%  
    filter(!( Code == "SQ" & str_detect(Value, "SEQUENCE" ))) %>%  
    group_by(Code, Gene) %>%  
    summarise(Value=str_c(Value, collapse = " ")) %>%  
    spread(Code, Value) %>%  
    select(Accessions=AC, Gene, Description=DE, GeneSymbol=GN, Sequence=SQ) %>%  
    mutate(Accessions=str_remove(Accessions, "; *$")) %>% 
    separate_rows(Accessions, sep="; ") %>%  
    mutate(GeneSymbol=str_remove_all(GeneSymbol, "^[[:alpha:]]+=|;.*")) %>%  
    mutate(Description=str_remove_all(Description, "RecName: |;$")) %>%  
    mutate(Description=str_remove_all(Description, "^Full=|;.*$")) %>%  
    mutate(Sequence=str_remove_all(Sequence, " "))
saveRDS(spTab, "annotation/SwissProtAnnotationTable.rds")
