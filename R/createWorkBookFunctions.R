# add worksheet with plots
addPlots <- function(wb, plotsObj){
    addWorksheet(wb, sheetName = "Coverage Plots", zoom = 75)
    showGridLines(wb, sheet = 1, showGridLines = FALSE)
    rowNum <- 1
    for(plotNum in seq_along(plotsObj)){
        print(plotsObj[[plotNum]][[1]])
        insertPlot(wb, 1, width = 22, height = 2, xy = NULL, startRow = rowNum,
               startCol = 1, fileType = "png", units = "in", dpi = 300)
        rowNum <- rowNum + 10
        print(plotsObj[[plotNum]][[2]])
        figHeight <- ceiling(width(plotsObj[[plotNum]][[3]]) / 100 ) * 0.6 + 0.2
        insertPlot(wb, 1, width = 22, height = figHeight, xy = NULL, startRow = rowNum,
               startCol = 1, fileType = "png", units = "in", dpi = 299)
        rowNum <- rowNum + (figHeight * 5)
    }
    wb
}
    
# header style for the tables
hs <- function(){
    createStyle(fontColour = "#ffffff", fgFill = "#4F80BD",
                  halign = "center", valign = "center",
                  border = "TopBottomLeftRight")
}

# Full data tables for each sample
addFullTables <- function(wb, protData){
    headStyle <- hs()
    for(sheetNum in seq_along(protData)){
        addWorksheet(wb, sampleTable$SampleName[sheetNum])
        showGridLines(wb, sheet = sheetNum, showGridLines = TRUE)
        writeData(wb, 
                  sheet = sheetNum + 1, 
                  x = protein_data[[sheetNum]], 
                  headerStyle = headStyle)
        setColWidths(wb, sheet = sheetNum + 1, cols = 1:3, widths = 16)
        setColWidths(wb, sheet = sheetNum + 1, cols = 4, widths = 100)
        setColWidths(wb, sheet = sheetNum + 1, cols = 5:6, widths = 18)
    }
    wb
}

# Combined data table - full table
addCombinedTable <- function(wb, combTab, sheetNum){
    headStyle <- hs()
    addWorksheet(wb, "Combined Data")
    writeData(wb, sheet = sheetNum, x = combTab, headerStyle = headStyle)
    setColWidths(wb, sheet = sheetNum, cols = 1:3, widths = 16) 
    setColWidths(wb, sheet = sheetNum, cols = 4, widths = 25) 
    setColWidths(wb, sheet = sheetNum, cols = 5:ncol(combTab), widths = "auto") 
    wb
}

# Filtered data table
addFilteredTable <- function(wb, filtTab, sheetNum){
    headStyle <- hs()
    addWorksheet(wb, "Filtered Data")
    writeData(wb, sheet = sheetNum, x = filtTab, headerStyle = headStyle)
    setColWidths(wb, sheet = sheetNum, cols = 1:3, widths = 16) 
    setColWidths(wb, sheet = sheetNum, cols = 4, widths = 25) 
    setColWidths(wb, sheet = sheetNum, cols = 5:ncol(filtTab), widths = "auto") 
    wb
}
    

makeWorkBook <- function(outputFile, plotsList, protData, combData, filtData){
                         
    # Create empty workbook
    wkbk <- createWorkbook(outputFile)
    # Generate Plots worksheet
    wkbk <- addPlots(wkbk, plotsList)
    # Add full per sample results tables
    wkbk <- addFullTables(wkbk, protData)
    # Add worksheet with combined table
    nSht <- length(protData) + 2
    wkbk <- addCombinedTable(wkbk, combData, nSht)
    # add worksheet with results filtgered by IgG
    nSht <- length(protData) + 3
    wkbk <- addFilteredTable(wkbk, filtData, nSht)
    
    # save workbook
    saveWorkbook(wkbk, outputFileName, overwrite = TRUE)
}
