library(shiny)

source("greg_power_script.r")
# fp005Text <- ""
# fp01Text <- ""
# powerData <- "uninitialized"
# nFalsePositives005 <- ""
# nFalsePositives01 <- ""

shinyServer(function(input, output, session) {

  


  powerData <- eventReactive(input$goButton, {
      output$fp005Text <- "Running simulation ..."
      powerData <- getHeatmapAndFDR(input$n)
      output$heatmap <- renderPlot({
        testPValues <- powerData$testPValues
        heatmapPalette <- colorRampPalette(c("deeppink3", "darkorchid1", "gray66"))(n = 299)
        col_breaks = c(seq(0,0.05,length=100),  # pink is significant at 0.05
          seq(0.05,0.1,length=100), #orchid is significant at 0.1
          seq(0.1,1,length=100)) #the rest is gray
        heatmap.2(testPValues,
                      cellnote = round(testPValues,digits=4),  # same data set for cell labels
                      main = "P Values", # heat map title
                      notecol="black",      # change font color of cell labels to black
                      density.info="none",  # turns off density plot inside color legend
                      trace="none",         # turns off trace lines inside the heat map
                     # margins =c(12,9),     # widens margins around plot
                      col=heatmapPalette,       # use on color palette defined earlier 
                      dendrogram='none',
                      breaks=col_breaks, 
                      Rowv=FALSE,
                      Colv=FALSE,
                      xlab="FOLD CHANGE",
                      ylab="BASE ABUNDANCE")
      })
      output$fp005Text <- "Number of false positives at 0.05 q-value with Benjamini-Hochberg multiple test correction:"
      output$nFalsePositives005 <- nrow(powerData$falsePositives005)
      output$fp01Text <- "Number of false positives at 0.1 q-value with Benjamini-Hochberg multiple test correction:"
      output$nFalsePositives01 <- nrow(powerData$falsePositives01)
      powerData
    })

  # output$fp005Text <- fp005Text
  # output$nFalsePositives005 <- nFalsePositives005

  # output$fp01Text <- fp01Text
  # output$nFalsePositives01 <- nFalsePositives01

})
