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
  titlePanel("Simulate power for a hypothetical count table experiment"),
  sidebarPanel(
        h4("Number of samples"),
        sliderInput("n", "Number of samples in control or experimental condition",
        2, 100, 5, step = 1)
    ),
  mainPanel(
    fluidRow(
      column(12,
        plotOutput("heatmap")
      )
    ),
    fluidRow(
      wellPanel(
          span("Number of false positives at 0.05 q-value with Benjamini-Hochberg multiple test correction:",
            textOutput("nFalsePositives005")
          )
          span("Number of false positives at 0.1 q-value with Benjamini-Hochberg multiple test correction:",
            textOutput("nFalsePositives01")
          )
        )
    )
  )
))
