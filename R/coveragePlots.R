#' Create fasta table
#' 
#' Creates a table with the fasta sequence for each bait protein for the report
#' @name createFastaTable
#' 
#' @import dplyr
#' @import purrr
#' @import readr
#' @importFrom Biostrings AAStringSet
createFastaTable <- function(sampleTab, annotTab){
    sampleTab %>% 
        select(BaitProteinUniProtID) %>% 
        distinct() %>% 
        filter(!is.na(BaitProteinUniProtID)) %>% 
        left_join(annotTab, by=c("BaitProteinUniProtID"="Accession")) %>% 
        mutate(ProteinSeq = map(Sequence, AAStringSet)) %>% 
        select(BaitProteinUniProtID, ProteinSeq)
}

#' Make summary coverage plot
#' 
#' This function returns takes the proteins amino acid sequence and the list of
#' overlapped amino acids (features) and generates a summary plot showing 
#' the overall coverage across the length of the protein.
#' @name makeCoveragePlots
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import purrr
#' @import stringr
covPlot1 <- function(Protein_seq, features, plotTitle){
    # get percent coverage
    protWidth <- str_length(as.character(Protein_seq))
    coverage <-  map2(features$start, features$end, seq) %>%   
        unlist() %>% 
        unique() %>% 
        length()
    Perct <- round(coverage / protWidth * 100, 2)
    SubTitle <- str_c("Number of Unique Peptides: ", 
                       nrow(features), 
                       "\n% Coverage: ", 
                       Perct)
    
    # set the tick positions for the plot
    nTicks <- min(c(7, ceiling(protWidth / 50) + 1))
    brkTicks <- round(seq(0, protWidth, length.out = nTicks), 0)
    
    features %>%
        distinct() %>%
        ggplot() +
        geom_rect(aes(xmin = start - 1, xmax = end, ymin = 0, ymax = 10), 
                  fill = "brown") +
        geom_rect(xmin = 0,
                  xmax = protWidth,
                  ymin = 0,
                  ymax = 10,
                  colour = "black",
                  fill = NA, 
                  size = 1) +
        labs(title = plotTitle, subtitle = SubTitle) +
        theme(
            axis.line=element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_rect(fill = "white"),
            panel.border = element_rect(colour = "black", fill=NA, size=3),
            plot.subtitle = element_text(hjust = 0.5),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_x_continuous(limits = c(0, protWidth), 
                           breaks = brkTicks,
                           expand = c(0,0)) +
        scale_y_continuous(limits = c(0, 10), 
                           breaks = c(0, 10), 
                           expand = c(0, 0))
}

#' Make detailed coverage plot
#' 
#' This function returns takes the proteins amino acid sequence and the list of
#' overlapped amino acids (features) and generates a detailed plot showing 
#' exactly which AAs are overlapped by peptides detected in the sample.
#' @name makeCoveragePlots
#' @import dplyr
#' @importFrom cowplot theme_nothing
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import purrr
#' @import stringr
covPlot2 <- function(Protein_seq, features){
    aaSeq <- str_split(as.character(Protein_seq), "")[[1]]
    len <- length(aaSeq)
    xIndex  <- rep(1:109, 100)
    xIndex <- xIndex[!(xIndex%%11==0)]
    xIndex <- xIndex[1:len]
    yIndex <- rep(1:100, each=100)[1:len]
    aaPosition <- 1:1e5
    aaPosition[!(aaPosition%%10%in%c(0, 1))] <- "..." 
    aaPosition <- aaPosition[1:len]
    aaPosition[len] <- len
    aaHits <- rep("", len)
    aaCovered <- map2(features$start, features$end, `:`) %>%   
        unlist()
    aaHits[aaCovered] <- "Covered"
    plotDat <- data.frame(Seq=aaSeq, 
                          Position=aaPosition,
                          xIndex=xIndex,
                          yIndex=yIndex,
                          highlight=aaHits)
    protWidth <- str_length(as.character(Protein_seq))
    yMax <- ceiling(protWidth / 100) + 1
    ggplot(plotDat) +
        geom_label(aes(x=xIndex, y=yIndex, fill=highlight, label=Seq), 
                   family = "mono") +
        geom_text(aes(x=xIndex, y=yIndex-0.4, label=Position)) +
        scale_fill_manual(values=c("#EEEDD6", "#A4F488")) +
        scale_y_reverse() +
        coord_cartesian(xlim = c(0.6, 109.4)) +
        theme_nothing() 
}

#' Create a feature table
#' 
#' Creates a table showing ranges of amino acids in the protein sequence that
#' are overlapped by the peptides detected in the sample
#' @name createFastaTable
#' @importFrom Biostrings vmatchPattern
#' @import dplyr
#' @importFrom magrittr %>%
getPosition <- function(peptideSeq, ProteinSeq) {
    vmatchPattern(peptideSeq, ProteinSeq) %>%
        as.data.frame() %>%
        select(start, end)
}

#' Make coverage plots for a single sample
#'
#' This function takes the bait protein sequences and peptide sequence data and
#' generates a table detailing which amino acids are overlaped by peptides
#' detected in the sample. It returns a list containing the summary coverage
#' map, the detailed coverage plot and the sequence for the bait protein.
#' @name coveragePlots
#' @import dplyr
#' @importFrom magrittr %>%
#' @import purrr
#' @import stringr
coveragePlots <- function(SampleName, BaitProteinName, BaitProteinUniProtID, 
                         PeptideData, ProteinSeq, ...) {

    features <- PeptideData %>%
        filter(`Master Protein Accessions` == BaitProteinUniProtID) %>%
        pull(Sequence) %>%
        toupper() %>%
        unique() %>%
        map_dfr(getPosition, ProteinSeq=ProteinSeq)
    
    title <- str_c(SampleName, ": ", BaitProteinName, " (", BaitProteinUniProtID, ")")
    cvp1 <- covPlot1(ProteinSeq, features, title)
    cvp2 <- covPlot2(ProteinSeq, features)
    list(CoverageMap=cvp1, 
         CoverageDetails=cvp2, 
         ProteinSeq=ProteinSeq)
}

#' Make coverage plots for each sample
#' 
#' This function retrieves the amino acid sequence for each bait protein, and 
#' then uses the information in the peptide tables to create two coverage plots
#' for each bait protein. The first shows the summary of the coverage as a bar,
#' the second shows the details of which amino acid sequences are covered by 
#' peptides detected.
#' @name makeCoveragePlots
#' @import dplyr
#' @import purrr
#' @importFrom magrittr %>%
makeCoveragePlots <- function(sampleTab, peptideDat, annotTab){
    fastaTab <- createFastaTable(sampleTab, annotTab)
    sampleTab %>%   
        mutate(PeptideData=peptideDat) %>% 
        inner_join(fastaTab, by="BaitProteinUniProtID") %>% 
        pmap(coveragePlots)
}
