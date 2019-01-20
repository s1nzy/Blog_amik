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
  titlePanel("PMX Solutions - Create a Shiny app for amikacin simulations"),
  
  # Sidebar 
  sidebarLayout(
    sidebarPanel(
      h4("Dosing details"),
      hr(),
      selectInput('admin', label='Route of administration', c("I.V. infusion","I.V. bolus")),
      numericInput("dose", label='Dose (mg/kg)', value = 12,min=0),
      numericInput("d2", label='Infusion duration (min)', value = 20,min=0),
      numericInput("ii1", label='Dosing interval (hours)', value = 24,min=0),
      
      
      br(),
      h4("Simulation details"),
      hr(),
      numericInput("nsim", label='Number of samples for simulation', value = 1000,min=0),
      numericInput("simtime", label='Simulation time', value = 125,min=0),
      numericInput("minprob", label='Minimal probability', value = 0.25,min=0,max=1),
      numericInput("maxprob", label='Maximal probability', value = 0.75,min=0,max=1),
      
      br(),
      h4("Patient covariates"),
      hr(),
      numericInput("bw", label='Birth weight (g)', value = 1730,min=0),
      numericInput("cw", label='Current weight (g)', value = 2800,min=0),
      numericInput("pna", label='Postnatal age (days)', value = 5,min=0),
      selectInput("nsaid", h4("Ibuprofen co-administration?"), 
                  choices = list("Yes" = 1, "No" = 0), selected = 0), 
      br(),
      h4("IIV"),
      hr(),
      numericInput("etav1", label='ETA - Central volume', value = 0.09,min=0),
      numericInput("etav2", label='ETA - Peripheral volume', value = 0,min=0),
      numericInput("etaq1", label='ETA - Inter-compartmental clearance', value = 0,min=0),
      numericInput("etacl", label='ETA - Clearance', value = 0,min=0),
      
      br(),
      h4("Residual variability"),
      hr(),
      numericInput("sigmaprop", label='Proportional', value = 0.0614,min=0),
      numericInput("sigmaadd", label='Additive', value = 0.267,min=0)
      
    ),

    mainPanel(
       plotOutput("simPlot"),
       downloadButton("downloadData", "Download")
       
    )
  )
))
