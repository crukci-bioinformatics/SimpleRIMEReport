### LOAD DATA FUNCTION #########################################################

# Columns to select from the Excel data tables
selectColumns <- function(x){
    select(x, 
           Accession, 
           Coverage, 
           `Unique Peptides`=`# Unique Peptides`)
}

# load protein data. Each table into a element of a list
loadProteinData <- function(sampleTab, dataDir){
    str_c(dataDir, "/", sampleTab$ProteinFile) %>% 
        map(read_tsv, col_types = cols()) %>% 
        map(selectColumns) %>% 
        set_names(sampleTab$SampleName)
}

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

checkBaitDetection <- function(sampleTab, dataDir){
    ProtData <- loadProteinData(sampleTab, dataDir)
    areGood <- sampleTab %>% 
        select(SampleName, BaitProteinName, BaitProteinUniProtID) %>% 
        mutate(ProteinDat=ProtData) %>%   
        pmap(checkProteins) %>%   
        unlist()
    sampleTab[areGood,]
}

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

getProteinTables <- function(sampleTab, dataDir, annot){
    loadProteinData(sampleTab, dataDir) %>%   
        map(left_join, annot, by=c("Accession"="Accessions")) %>% 
        map(select, Accession, Gene, `Gene Symbol`=GeneSymbol, 
            Description, everything())
}

# load peptide data. Each table into a element of a list
getPeptideTables <- function(sampleTab, dataDir){
    str_c(dataDir, "/", sampleTab$PeptideFile) %>% 
        map(read_tsv, col_types = cols()) %>% 
        set_names(sampleTab$SampleName)
}

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
    

#    non_specific <- combTab %>% 
#        select(Accession, contains("IGG")) %>% 
#        gather("ColNam", "ANumber", -Accession) %>% 
#        filter(!is.na(ANumber)) %>% 
#        pull(Accession) %>% 
#        unique()
#        
#    combTab %>% 
#        filter(!Accession%in%non_specific) %>% 
#        select(-contains("IGG"))
#}
################################################################################
