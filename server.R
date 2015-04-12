#library(ggvis)
#library(dplyr)
library(shiny)
library(RSQLite)
library(RSQLite.extfuns)

# Set up handles to database tables on app start
db <- src_sqlite("sparsityData")

sparsity <- list()
for (i in 0:10) {
  sparsity[[i+1]] <- tbl(db,paste("sparsity",i*10,sep=""))
}

shinyServer(function(input, output, session) {

  plotdata <- reactive({
    as.data.frame(sparsity[[input$reviews/10 + 1]])
  })

  output$plot1 <- renderPlot({
    plotdata <- plotdata()
    stripchart(plotdata[,c(2,10,8,3,11,9)], main="% Var explained", method="jitter", vertical=T, group.names=c("P", "P.r", "P.D", "clr", "clr.r", "clr.D"), pch=19, col=rgb(0,0,0,0.3),ylim=c(0,1))
  })
  output$plot2 <- renderPlot({
    plotdata <- plotdata()
    stripchart(plotdata[,c(4,12,6,5,13,7)], main="Corr with first PC", method="jitter", vertical=T, group.names=c("P", "P.r", "P.D", "clr", "clr.r", "clr.D"), pch=19, col=rgb(0,0,0,0.3))
  })

  output$sparsityFilter <- renderText({ input$reviews })
})
