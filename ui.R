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
        2, 100, 5, step = 1),
        br(),
        p("Click the button to perform the power simulation. This may take a long time depending on how large your data set is."),
        actionButton("goButton", "Go!")
    ),
  mainPanel(
    fluidRow(
      column(12,
        plotOutput("heatmap")
      )
    ),
    fluidRow(
      wellPanel(
          span(textOutput("fp005Text"),
            textOutput("nFalsePositives005")
          ),
          span(textOutput("fp01Text"),
            textOutput("nFalsePositives01")
          )
        )
    )
  )
))
