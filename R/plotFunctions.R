#### PLOTTING FUNCTION #########################################################

## match sequence with protein sequence and return co-ordinates
getPosition <- function(peptideSeq, ProteinSeq) {
    vmatchPattern(peptideSeq, ProteinSeq) %>%
        as.data.frame() %>%
        dplyr::select(start, end) %>%
        return()
}

covPlot1 <- function(Protein_seq, features, plotTitle){
    # get percent coverage
    protWidth <- width(Protein_seq)
    coverage <- GRanges("feature", IRanges(features$start, features$end)) %>%
        GenomicRanges::reduce() %>%
        width() %>%
        sum()
    Perct <- round(coverage / protWidth * 100, 2)
    SubTitle <- paste0("Number of Unique Peptides: ", 
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
                           breaks = brkTicks) +
                           #expand = c(0,0)) +
        scale_y_continuous(limits = c(0, 10), 
                           breaks = c(0, 10), 
                           expand = c(0, 0))
}

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
    yMax <- ceiling(width(Protein_seq) / 100) + 1
    ggplot(plotDat) +
        geom_label(aes(x=xIndex, y=yIndex, fill=highlight, label=Seq), 
                   family = "mono") +
        geom_text(aes(x=xIndex, y=yIndex-0.4, label=Position)) +
        scale_fill_manual(values=c("#EEEDD6", "#A4F488")) +
        scale_y_reverse() +
        coord_cartesian(xlim = c(0.6, 109.4)) +
        #coord_cartesian(xlim = c(0.6, 109.4), ylim = c(0, yMax)) +
        theme_nothing() 
}

# Coverage plot
coveragePlot <- function(SampleName, BaitProteinName, BaitProteinUniProtID, 
                         PeptideData, ProteinSeq, species_id) {
    message("Generate plots for ", SampleName)

    features <- PeptideData %>%
        filter(`Protein Group Accessions` == BaitProteinUniProtID) %>%
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

# make Plots
makeCoveragePlots <- function(sampleTab, peptideDat, fastaTab){
    sampleTab %>%   
        mutate(PeptideData=peptideDat) %>% 
        inner_join(fastaTab) %>% 
        select(SampleName, BaitProteinName, 
               BaitProteinUniProtID, PeptideData, ProteinSeq) %>% 
        pmap(coveragePlot, species_id = species_id)
}
