#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

###############################################
############ Load libraries
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  
  

  
  
  
  sim_dataframe <- reactive({
    
    ###############################################
    ### Dosing
    if(input$admin =='Depot'){
      oral <- T
    }else{
      oral <- F
    }
    
    dose <- input$dose
    
    
    ###############################################
    ## Parameter estimates 
    ka  <- input$ka # Absorption rate constant
    cl  <- input$cl # Clearance
    vd  <- input$v1 # Central distribution volume
    vd2 <- input$v2 # Peripheral distribution volume
    q1  <- input$q1 # Inter-compartmental clearance (set to 0 for 1-cmt model)
    
    
    
    
    ###############################################
    ## Insert etas and sigmas
    
    etaka  <- input$etaka
    etacl  <- input$etacl
    etavd  <- input$v1
    etavd2 <- input$v2
    etaq1  <- input$q1
    
    
    
    
    #################################################
    ## Insert residual error
    
    sigmaprop <- input$sigmaprop # Proportional error
    sigmaadd  <- input$sigmaadd # Additive error
    
    
    
    
    ###############################3
    ## Simulation info
    
    nsamples <- input$nsim ### Number of simulated individuals
    sim_time <- input$simtime ## Time of simulation
    
    
    

    ###################################################
    ############# Set Dosing objects
    
    if(oral){
      ## Oral dose
      Administration <-  as.data.frame(ev(ID=1:nsamples,ii=24, cmt=1, addl=9999, amt=dose, rate = 0,time=0)) 
    }else{
      ## IV BOLUS
      Administration <-  as.data.frame(ev(ID=1:nsamples,ii=24, cmt=2, addl=9999, amt=dose, rate = 0,time=0)) 
    }
    
    
    
    ## Sort by ID and time
    data <- Administration[order(Administration$ID,Administration$time),]
    
    
    
    ######################################################################
    ### Load in model for mrgsolve
    mod <- mread_cache("popPK")
    
    
    
    ## Specify the omegas and sigmas matrices
    omega <- cmat(etaka,
                  0, etacl,
                  0,0,etavd,
                  0,0,0,etavd2,
                  0,0,0,0,etaq1)
    
    
    sigma <- cmat(sigmaprop,
                  0,sigmaadd)
    
    
    
    ## Set parameters in dataset
    data$TVKA <- ka
    data$TVCL <- cl
    data$TVVC <- vd
    data$TVVP1<- vd2
    data$TVQ1 <- q1
    
    
    
    
    #################################################################
    ###### Perform simulation
    out <- mod %>%
      data_set(data) %>%
      omat(omega) %>%
      smat(sigma) %>%
      mrgsim(end=sim_time,delta=sim_time/100, obsonly=TRUE) # Simulate 100 observations, independent of the total simulated time
    
    
    ### Output of simulation to dataframe
    df <- as.data.frame(out)
    
    
    return(df)
  })
  

  
  
  
sum_stat <- reactive({
    
  #################################################################  
  ## Calculate summary statistics
  
  
  # set probabilities of ribbon in figure
  minprobs=input$minprob
  maxprobs=input$maxprob
  
  
  sum_stat <- sim_dataframe() %>%
    group_by(time) %>%
    summarise(Median_C=median(DV),
              Low_percentile=quantile(DV,probs=minprobs),
              High_percentile=quantile(DV,probs=maxprobs)
    )  
  return(sum_stat)
  
})
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  output$simPlot <- renderPlot({
    
    ggplot(sum_stat(), aes(x=time,y=Median_C)) +
      
      ## Add ribbon for variability
      geom_ribbon(aes(ymin=Low_percentile, ymax=High_percentile, x=time), alpha = 0.15, linetype=0)+
      
      ## Add median line
      geom_line(size=2) +
      
      # scale_y_log10()+
      
      # Set axis and theme
      ylab(paste("Concentration",sep=""))+
      xlab("Time after dose (h)")+
      theme_bw()+
      
      # Remove legend
      theme(legend.position="none")
    
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      "Simulated_dataset.csv"
    },
    content = function(file) {
      write.csv(sim_dataframe(), file, row.names = FALSE)
    }
  )
  
})