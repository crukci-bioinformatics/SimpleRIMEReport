#### CREATING THE INITIAL SAMPLE TABLE #########################################

# split the file name and extract the sample name and protein
filenametotable <- function(fileName){
    dat <- str_split(fileName, "_")[[1]]
    len <- length(dat)
    tibble(ProteinFile=fileName,
           PeptideFile=str_replace(ProteinFile, "proteingroups", "psms"), 
           SampleName=str_c(dat[6:(len-2)], collapse="_"),
           BaitProteinName=dat[len-1])
}
# create the sample table from files in the data directory
createSampleTable <- function(datDir){
    list.files(datDir, pattern="proteingroups.txt$") %>%
        map_df(filenametotable) %>% 
        mutate(BaitProteinUniProtID="") %>% 
        write_csv(str_c(datDir, "/", basename(datDir), "_sample_details.csv"))
}
