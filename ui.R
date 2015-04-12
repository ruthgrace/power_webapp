library(ggvis)
library(shiny)

# For dropdown menu
actionLink <- function(inputId, ...) {
  tags$a(href='javascript:void',
         id=inputId,
         class='action-button',
         ...)
}

shinyUI(fluidPage(
  titlePanel("Variance and Read Count"),
  sidebarPanel(
        h4("Sparsity Filter"),
        sliderInput("reviews", "Minimum number of reads per OTU",
        0, 100, 50, step = 10)
    ),
  mainPanel(
    fluidRow(
      column(6,
        plotOutput("plot1")
      ),
      column(6,
        plotOutput("plot2")
      )
    ),
    fluidRow(
      wellPanel(
          span("Minimum counts per OTU across all samples:",
            textOutput("sparsityFilter")
          )
        )
    )
  )
))
