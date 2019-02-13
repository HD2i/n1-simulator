#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Author: Michael Jones
# Organization: HD2i @ INGH @ MSSM

library(shiny)
library(ggplot2)
source('../n1-simulator.R')
source('utils.R')

# Define UI for application
ui <- function(request) {
  fixedPage(
    includeCSS("custom.css"),
    
    title = "N-of-1 Study Simulator",
    titlePanel("N-of-1 Study Simulator"),
    hr(),
    
    fixedRow(id="datavisRow",
      column(12,
        tabsetPanel(
          tabPanel("Plot", 
            plotOutput('outcomePlot')
          ),
          tabPanel("Fit", 
            h4("Model Fit 1: Treatment Only"),
            verbatimTextOutput("model1fit"),
            hr(),
            h4("Model Fit 2: Treatment + Block"),
            verbatimTextOutput("model2fit"),
            hr(),
            h4("Model Fit 3: Treatment + Block + Treatment*Time"),
            verbatimTextOutput("model3fit"),
            hr(),
            h4("Model Fit 4: Treatment + Block Mixed Effect"),
            verbatimTextOutput("model4fit")
          )
        )
      )
    ),
    hr(),
    fixedRow(id="inputsRow",
      column(4,
        h3("Treatments"),
        sliderInput("n_treatments", "Number of treatments:", value=2, min=1, max=5),
        hr(),
        h4("Treatment A"),
        numericInput("effect_1", "Effect Size:", value=-40),
        numericInput("tc_in_1", "Run In:", min=1, value=3),
        numericInput("tc_out_1", "Wash Out:", min=1, value=2),
        conditionalPanel(
         condition = "input.n_treatments >= 2",
         hr(),
         h4("Treatment B"),
         numericInput("effect_2", "Treatment B Effect Size:", value=-30),
         numericInput("tc_in_2", "Treatment B Run In:", min=1, value=2),
         numericInput("tc_out_2", "Treatment B Wash Out:", min=1, value=3)
        ),
        conditionalPanel(
         condition = "input.n_treatments >= 3",
         hr(),
         h4("Treatment C"),
         numericInput("effect_3", "Treatment C Effect Size:", value=-40),
         numericInput("tc_in_3", "Treatment C Run In:", min=1, value=3),
         numericInput("tc_out_3", "Treatment C Wash Out:", min=1, value=2)
        ),
        conditionalPanel(
         condition = "input.n_treatments >= 4",
         hr(),
         h4("Treatment D"),
         numericInput("effect_4", "Treatment D Effect Size:", value=-40),
         numericInput("tc_in_4", "Treatment D Run In:", min=1, value=3),
         numericInput("tc_out_4", "Treatment D Wash Out:", min=1, value=2)
        ),
        conditionalPanel(
         condition = "input.n_treatments >= 5",
         hr(),
         h4("Treatment E"),
         numericInput("effect_5", "Treatment E Effect Size:", value=-40),
         numericInput("tc_in_5", "Treatment E Run In:", min=1, value=3),
         numericInput("tc_out_5", "Treatment E Wash Out:", min=1, value=2)
        )
      ),
      column(4,
        h3("Study Design"),
        sliderInput("n_blocks", "Number of treatment blocks:", value=2, min=1, max=10),
        uiOutput("treatmentScheduleOptions"),
        numericInput("treatment_period", "Treatment Period:", value=30),
        hr(),
        h3("Statistical Parameters"),
        numericInput("sd_baseline", "S.D. Baseline Drift", min=0, value=0.4, step=0.1),
        numericInput("sd_outcome", "S.D. Process Noise", min=0, value=0.6, step=0.1),
        numericInput("sd_obs", "S.D. Observation Noise", min=0, value=4, step = 0.1)
      ),
      column(4,
        h3("Randomization"),
        numericInput("random_seed", "Random Seed:", 1),
        actionButton("random_seed_generate", "Generate"),
        hr(),
        actionButton("reset_button", "Reset inputs to default"),
        br(),
        br(),
        bookmarkButton()
      )
    ),
    hr(),
    fixedRow(id="orgRow",
      column(4,
        a(img(src='hd2ilogo.svg', height=40), href="http://hd2i.org")
      ),
      column(4,
        a(img(src='inghlogo.svg', height=60), href="http://www.nextgenhealthcare.org")
      ),
      column(4,
        a("Source available on Github", href="https://github.com/HD2i/n1-simulator")
      )
    )
  )
}

# Define server logic
server <- function(input, output, session) {
  ## Notification
  notification_id <- NULL
  
  observe({
    input$n_treatments
    input$n_blocks
    input$effect_1
    input$effect_2
    input$effect_3
    input$effect_4
    input$effect_5
    input$tc_in_1
    input$tc_in_2
    input$tc_in_3
    input$tc_in_4
    input$tc_in_5
    input$tc_out_1
    input$tc_out_2
    input$tc_out_3
    input$tc_out_4
    input$tc_out_5
    input$treatment_period
    input$sd_baseline
    input$sd_outcome
    input$sd_obs
    input$random_seed
    notification_id <<- showNotification("Working...", duration = NULL, closeButton=FALSE, type="message")
  })
  
  observeEvent(n1sim_result(), {
    if (!is.null(notification_id))
      removeNotification(notification_id)
    notification_id <<- NULL
  })
  
  ## Reactive elements
  effect_sizes <- reactive({
    req(input$treatment_schedule)
    input$treatment_schedule
    
    result <- vector()
    for(i in 1:isolate(input$n_treatments)) {
      result <- c(result, eval(parse(text=paste0("input$effect_",i))))
    }
    return(result)
  })
  
  tc_ins <- reactive({
    req(input$treatment_schedule)
    input$treatment_schedule
    
    result <- vector()
    for(i in 1:isolate(input$n_treatments)) {
      result <- c(result, eval(parse(text=paste0("input$tc_in_",i))))
    }
    return(result)
  })
  
  tc_outs <- reactive({
    req(input$treatment_schedule)
    input$treatment_schedule
    
    result <- vector()
    for(i in 1:isolate(input$n_treatments)) {
      result <- c(result, eval(parse(text=paste0("input$tc_out_",i))))
    }
    return(result)
  })
  
  n1sim_result <- reactive({
    n1_simulate(
      n_treatments = isolate(input$n_treatments), 
      n_blocks = isolate(input$n_blocks),
      baseline_initial = 160,
      effect_size = effect_sizes(),
      tc_in = tc_ins(),
      tc_out = tc_outs(),
      tc_outcome = 1,
      sd_baseline = input$sd_baseline,
      sd_outcome = input$sd_outcome,
      sd_obs = input$sd_obs,
      treatment_period = input$treatment_period,
      sampling_timestep = 1.0,
      noise_timestep = 0.01,
      treatment_order = treatment_order_str_to_mat(isolate(input$n_blocks),input$treatment_schedule),
      random_seed = input$random_seed,
      return_data_frame = TRUE
    )
  })
  
  treatment_options <- reactive({
    treatment_schedule_options_to_strvec(generate_treatment_schedule_options(input$n_treatments, input$n_blocks))
  })
  
  ## Render elements
  output$treatmentScheduleOptions <- renderUI({
    selectInput("treatment_schedule", "Treatment Schedule:", treatment_options())
  })
  
  output$outcomePlot <- renderPlot({
    # prevent plotting before treatment_schedule has been set
    req(input$treatment_schedule)
    
    timeseries <- n1sim_result()$timeseries
    plot <- ggplot(data = timeseries, aes(x = t)) 
    plot <- plot + geom_line(aes(y = outcome, alpha="underlying"), colour="#555555", size=0.75)
    plot <- plot + geom_point(aes(y = outcome_obs, alpha="observed", colour=factor(treatment), shape=factor(treatment), fill=factor(treatment)))
    
    isolate({
      grid_x_major = seq(0, input$treatment_period*input$n_blocks*input$n_treatments, input$treatment_period)
      grid_x_minor = seq(input$treatment_period/2, input$treatment_period*(input$n_blocks*input$n_treatments-1/2), input$treatment_period)
    })
    plot <- plot + scale_x_continuous(minor_breaks = grid_x_minor, breaks = grid_x_major)
    plot <- plot + ggtitle("Underlying and observed outcome") + theme(plot.title = element_text(hjust=0.5, face="bold"))
    plot <- plot + xlab("Day in study") + ylab("Outcome")
    plot <- plot + theme(text = element_text(size=14))
    plot <- plot + scale_colour_manual(name="Treatment",
                                       values=c("1"="#009E73","2"="#D55E00","3"="#0072B2","4"="#CC79A7","5"="#E69F00"))
    plot <- plot + scale_fill_manual(name="Treatment",
                                     values=c("1"="#009E73","2"="#D55E00","3"="#0072B2","4"="#CC79A7","5"="#E69F00"))
    plot <- plot + scale_shape_manual(name="Treatment",
                                      values=c("1"=21,"2"=22,"3"=23,"4"=24,"5"=25))
    plot <- plot + scale_alpha_manual(name="Outcome",
                                      values=c("observed"=1,"underlying"=1),
                                      guide=guide_legend(override.aes = list(linetype=c("blank","solid"),shape=c(16,NA))))
    return(plot)
  })
  
  output$model1fit <- renderPrint({
    n1_fit_model1(n_treatments = isolate(input$n_treatments), n1sim_result()$timeseries)
  })
  output$model2fit <- renderPrint({
    n1_fit_model2(n_treatments = isolate(input$n_treatments), n1sim_result()$timeseries)
  })
  output$model3fit <- renderPrint({
    n1_fit_model3(n_treatments = isolate(input$n_treatments), n1sim_result()$timeseries)
  })
  output$model4fit <- renderPrint({
    n1_fit_model4(n_treatments = isolate(input$n_treatments), n1sim_result()$timeseries)
  })
  
  ## Observer elements
  observeEvent(input$random_seed_generate, {
    updateNumericInput(session, "random_seed", value = sample(2^31 - 1, 1))
  })
  
  observeEvent(input$reset_button, {
    updateNumericInput(session, "n_treatments", value = 2)
    updateNumericInput(session, "n_blocks", value = 2)
    
    updateNumericInput(session, "effect_1", value = -40)
    updateNumericInput(session, "tc_in_1", value = 3)
    updateNumericInput(session, "tc_out_1", value = 2)
    updateNumericInput(session, "effect_2", value = -30)
    updateNumericInput(session, "tc_in_2", value = 2)
    updateNumericInput(session, "tc_out_2", value = 3)
    updateNumericInput(session, "effect_3", value = -40)
    updateNumericInput(session, "tc_in_3", value = 3)
    updateNumericInput(session, "tc_out_3", value = 2)
    updateNumericInput(session, "effect_4", value = -40)
    updateNumericInput(session, "tc_in_4", value = 3)
    updateNumericInput(session, "tc_out_4", value = 2)
    updateNumericInput(session, "effect_5", value = -40)
    updateNumericInput(session, "tc_in_5", value = 3)
    updateNumericInput(session, "tc_out_5", value = 2)

    updateSliderInput(session, "treatment_period", value = 30)
    
    updateNumericInput(session, "sd_baseline", value = 0.4)
    updateNumericInput(session, "sd_outcome", value = 0.6)
    updateNumericInput(session, "sd_obs", value = 4)
    
    updateNumericInput(session, "random_seed", value = 1)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server, enableBookmarking="url")

