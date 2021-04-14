library(gridExtra)

numberOfRows <- function(plotObj){
    ceiling(width(plotObj$ProteinSeq) / 100)
}

makeDualPlot <- function(plotObj, maxRow){
    thisRow <- numberOfRows(plotObj)
    figHeight2a <- thisRow * 0.6 + 0.2
    figHeight2b <- (maxRow * 0.6 + 0.2) - figHeight2a
    emptyPlot <- ggplot() + theme(axis.line=element_blank())
    grid.arrange(plotObj[[1]], 
                 plotObj[[2]], 
                 emptyPlot, 
                 nrow = 3, 
                 heights = c(1.8, figHeight2a, figHeight2b))
}

maxRow <- map_dbl(coveragePlots, numberOfRows) %>%  
    max()

pdf("TestData.pdf", width = 25, height=17)
coveragePlots %>%  
    walk(makeDualPlot, maxRow=maxRow)
dev.off()
    

