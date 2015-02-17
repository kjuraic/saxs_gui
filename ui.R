
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("SAXS Analysis"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      fileInput('saxsFile', 'Choose 1D SAXS file',
                accept=c('*.dat')),
      tags$hr(),
      sliderInput("rangePlotGun", "Guinier plot range:", min = 1, max = 1000, value = c(1,1000)),
      sliderInput("rangeGun", "Guinier fit range:", min = 1,max = 1000,value = c(1,1000)),
      sliderInput("rangePlotPor", "Porod plot range:", min = 1, max = 1000, value = c(1,1000)),
      sliderInput("rangePor", "Porod fit range:", min = 1, max = 1000, value = c(1,1000)),
      sliderInput("rangePlotInv", "Invariant plot range:", min = 1, max = 1000, value = c(1,1000)),
      sliderInput("rangeInv", "Invariant fit range:", min = 1, max = 1000, value = c(1,1000))
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs", 
                  tabPanel("SAXS plot", 
                           plotOutput("saxsPlot")),
                  tabPanel("Guinier"  , plotOutput("plotGun" ),
                                        tableOutput("summaryGun")
                  ),
                  tabPanel("Porod", plotOutput("plotPor"),
                                        tableOutput("summaryPor")), 
                  tabPanel("Invariant", plotOutput("plotInv")),
                  tabPanel("Summary", textOutput("SAXS summary"))
      )
    )
  )
))