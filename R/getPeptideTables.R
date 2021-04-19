
# load peptide data. Each table into a element of a list
getPeptideTables <- function(sampleTab, dataDir){
    str_c(dataDir, "/", sampleTab$PeptideFile) %>% 
        map(read_tsv, col_types = cols()) %>% 
        set_names(sampleTab$SampleName)
}