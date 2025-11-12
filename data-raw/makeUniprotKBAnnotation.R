# The purpose of the script is to create an annotation object from UniProt.
# This object is saved as an internal data source (R/sysdata.rda) and is used by
# the package annotate the data set. This saves having to use Uniprot.ws, which
# can be quite slow. It unforunately does create a very large data object.
# First you need to download SwissProt from:
# 
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
#
# This file is > 500Mb, so you should then filter it before running this script:
#    zcat uniprot_sprot.dat.gz |
#      grep -E "^ID|^AC|^OS|^GN|^DE|^SQ|^[[:blank:]]"  |
#      gzip -c - > uniprot_sprot.trim.dat.gz

library(tidyverse)

keepSpecies <- str_c("Homo sapiens",
                     "|Mus musculus",
                     "|Rattus",
                     "|Drosophila",
                     "|Caenorhabditis elegans")

spTab <- read_tsv("~/Downloads/uniprot_sprot.trim.dat.gz",
                  col_names = FALSE,
                  trim_ws = FALSE) %>%
    separate(X1, c("Code", "Value"), extra = "merge")

annot <- spTab %>%
    filter(Code %in% c("ID", "AC", "GN", "DE", "OS", "SQ", "")) %>%
    filter(Value != "") %>%
    mutate(Gene = ifelse(Code == "ID", Value, NA)) %>%
    mutate(Gene = str_remove(Gene, " .*")) %>%
    fill(Gene, .direction = "down") %>%
    filter(!(Code == "DE" & str_detect(Value, "RecName", negate = TRUE)))  %>%
    mutate(Code = str_replace(Code, "^$", "SQ")) %>%
    filter(!(Code == "SQ" & str_detect(Value, "SEQUENCE"))) %>%
    summarise(Value = str_c(Value, collapse = " "), .by = c("Code", "Gene")) %>%
    pivot_wider(names_from = Code, values_from = Value) %>%
    select(Accession = AC,
           Gene,
           Description = DE,
           GeneSymbol = GN,
           Sequence = SQ,
           Species = OS) %>%
    filter(str_detect(Species, keepSpecies)) %>%
    select(-Species) %>%
    mutate(Accession = str_remove(Accession, "; *$")) %>%
    separate_rows(Accession, sep = "; ") %>%
    mutate(GeneSymbol = str_remove_all(GeneSymbol,
                                       "^[[:alpha:]]+=|;.*| \\{ECO.*")) %>%
    mutate(Description = str_remove_all(Description, "RecName: |;$")) %>%
    mutate(Description = str_remove_all(Description, "^Full=|;.*$")) %>%
    mutate(Sequence = str_remove_all(Sequence, " "))

save(annot, file = "data/uniprot_annotation.rda", version = 3)
