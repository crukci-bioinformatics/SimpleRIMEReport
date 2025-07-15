#' Add worksheet with plots
#'
#' This function adds the coverage plots to a new worksheet in the workbook. The
#' summary coverage map is always the same height, but the plot height for the
#' detailed coverage map needs to be calculated from the amino acid sequence
#' length for the bait protein.
#' @name addPlots
#' @import openxlsx
addPlots <- function(wb, plotsObj){
    addWorksheet(wb, sheetName = "Coverage Plots", zoom = 75)
    showGridLines(wb, sheet = 1, showGridLines = FALSE)
    rowNum <- 1
    for(plotNum in seq_along(plotsObj)){
        # summary coverage map
        print(plotsObj[[plotNum]]$CoverageMap)
        insertPlot(wb,
                   sheet = 1,
                   width = 22,
                   height = 2,
                   xy = NULL,
                   startRow = rowNum,
                   startCol = 1,
                   fileType = "png",
                   units = "in",
                   dpi = 300)
        rowNum <- rowNum + 10

        # detailed coverage map
        print(plotsObj[[plotNum]]$CoverageDetails)
        protWidth <- str_length(as.character(plotsObj[[plotNum]]$ProteinSeq))
        figHeight <- ceiling(protWidth / 100 ) * 0.6 + 0.2
        insertPlot(wb,
                   sheet = 1,
                   width = 22,
                   height = figHeight,
                   xy = NULL,
                   startRow = rowNum,
                   startCol = 1,
                   fileType = "png",
                   units = "in",
                   dpi = 299)
        rowNum <- rowNum + (figHeight * 5)
    }
    wb
}

#' Header style for the tables
#' @name hs
#' @import openxlsx
hs <- function(){
    createStyle(fontColour = "#ffffff",
                fgFill = "#4F80BD",
                halign = "center",
                valign = "center",
                border = "TopBottomLeftRight")
}

#' Add full data tables for each sample
#'
#' This function adds one worksheet for each sample, with the complete protein
#' details
#' @name addFullTables
#' @import openxlsx
addFullTables <- function(wb, protData, nSht){
    headStyle <- hs()
    for(sheetNum in seq_along(protData)){
        addWorksheet(wb, names(protData)[sheetNum])
        showGridLines(wb, sheet = sheetNum + nSht, showGridLines = TRUE)
        writeData(wb,
                  sheet = sheetNum + nSht,
                  x = protData[[sheetNum]],
                  headerStyle = headStyle)
        setColWidths(wb, sheet = sheetNum + nSht, cols = 1:3, widths = 16)
        setColWidths(wb, sheet = sheetNum + nSht, cols = 4, widths = 100)
        setColWidths(wb, sheet = sheetNum + nSht, cols = 5:6, widths = 18)
    }
    wb
}

#' Add complete combined data table
#'
#' This function adds a worksheet with the combined table showing protein
#' details for all samples
#' @name addCombinedTable
#' @import openxlsx
addCombinedTable <- function(wb, combTab, sheetNum){
    headStyle <- hs()
    addWorksheet(wb, "Combined Data")
    writeData(wb, sheet = sheetNum, x = combTab, headerStyle = headStyle)
    setColWidths(wb, sheet = sheetNum, cols = 1:3, widths = 16)
    setColWidths(wb, sheet = sheetNum, cols = 4, widths = 25)
    setColWidths(wb, sheet = sheetNum, cols = 5:ncol(combTab), widths = "auto")
    wb
}

#' Add complete filtered data table
#'
#' This function adds a worksheet with the table showing protein details for all
#' samples but with the non-specific binding proteins filtered out
#' @name addFilteredTable
#' @import openxlsx
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

#' Create the workbook
#'
#' This function creates a new work book and then adds the plots, individual
#' data tables for each sample, the complete combined data and the filtered
#' combinded data
#' @name makeWorkBook
#' @import openxlsx
makeWorkBook <- function(outputFile, plotsList, protData, combData, filtData){

    # Create empty workbook
    wkbk <- createWorkbook(outputFile)
    nSht <- 0
    # Generate Plots worksheet
    if(length(plotsList) > 0) {
        wkbk <- addPlots(wkbk, plotsList)
        nSht <- nSht + 1
    } else {
        message("No coverage plots were generated, ",
                "skipping the plots worksheet.")
    }
    # Add full per sample results tables
    wkbk <- addFullTables(wkbk, protData, nSht)
    nSht <- length(protData) + nSht
    # Add worksheet with combined table
    wkbk <- addCombinedTable(wkbk, combData, nSht + 1)
    nSht <- nSht + 1
    # add worksheet with results filtgered by IgG
    wkbk <- addFilteredTable(wkbk, filtData, nSht + 1)

    # save workbook
    saveWorkbook(wkbk, outputFile, overwrite = TRUE)
}
