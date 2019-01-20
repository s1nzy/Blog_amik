#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("PMX Solutions - Create a Shiny app for PK simulations"),
  
  # Sidebar 
  sidebarLayout(
    sidebarPanel(
      h4("Dosing details"),
      hr(),
      selectInput('admin', label='Route of administration', c("Depot","I.V. bolus")),
      numericInput("dose", label='Dose', value = 100,min=0),
      br(),
      h4("Simulation details"),
      hr(),
      numericInput("nsim", label='Number of samples for simulation', value = 1000,min=0),
      numericInput("simtime", label='Simulation time', value = 12,min=0),
      numericInput("minprob", label='Minimal probability', value = 0.25,min=0,max=1),
      numericInput("maxprob", label='Maximal probability', value = 0.75,min=0,max=1),
      
      br(),
      h4("Population parameters"),
      hr(),
      numericInput("ka", label='Absorption rate constant', value = 10,min=0),
      numericInput("v1", label='Central volume', value = 10,min=0),
      numericInput("v2", label='Peripheral volume', value = 10,min=0),
      numericInput("q1", label='Inter-compartmental clearance', value = 10,min=0),
      numericInput("cl", label='Clearance', value = 10,min=0),
      
      br(),
      h4("IIV"),
      hr(),
      numericInput("etaka", label='ETA - Absorption rate constant', value = 10,min=0),
      numericInput("etav1", label='ETA - Central volume', value = 10,min=0),
      numericInput("etav2", label='ETA - Peripheral volume', value = 10,min=0),
      numericInput("etaq1", label='ETA - Inter-compartmental clearance', value = 10,min=0),
      numericInput("etacl", label='ETA - Clearance', value = 10,min=0),
      
      br(),
      h4("Residual variability"),
      hr(),
      numericInput("sigmaprop", label='Proportional', value = 0.01,min=0),
      numericInput("sigmaadd", label='Additive', value = 0.1,min=0)
      
    ),

    mainPanel(
       plotOutput("simPlot"),
       downloadButton("downloadData", "Download")
       
    )
  )
))
