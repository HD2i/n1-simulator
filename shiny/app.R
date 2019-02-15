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
library(shinyBS)
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
        numericInput("effect_1", "Effect size:", value=-40),
        numericInput("tc_in_1", HTML(paste0("Run-in ", tags$span("(days)"), ":")), min=1, value=3),
        numericInput("tc_out_1", HTML(paste0("Wash-out ", tags$span("(days)"), ":")), min=1, value=2),
        bsTooltip("effect_1", "The effect size is the asymptotic magnitude change from baseline that a person would experience in the long run, if they continued this treatment",
                  placement="right", trigger="hover", options=list(container="body")),
        bsTooltip("tc_out_1", "The \"wash-out\" (or \"carryover\") is a measure of how long it takes for a treatment to stop working once it is discontinued",
                  placement="right", trigger="hover", options=list(container="body")),
        bsTooltip("tc_in_1", "The \"run-in\" (or \"wash-in\") is a measure of how long it takes for a treatment to ramp up to its full effect",
                  placement="right", trigger="hover", options=list(container="body")),
        conditionalPanel(
         condition = "input.n_treatments >= 2",
         hr(),
         h4("Treatment B"),
         numericInput("effect_2", "Effect size:", value=-30),
         numericInput("tc_in_2", HTML(paste0("Run-in ", tags$span("(days)"), ":")), min=1, value=2),
         numericInput("tc_out_2", HTML(paste0("Wash-out ", tags$span("(days)"), ":")), min=1, value=3),
         bsTooltip("effect_2", "The effect size is the asymptotic magnitude change from baseline that a person would experience in the long run, if they continued this treatment",
                   placement="right", trigger="hover", options=list(container="body")),
         bsTooltip("tc_out_2", "The \"wash-out\" (or \"carryover\") is a measure of how long it takes for a treatment to stop working once it is discontinued",
                   placement="right", trigger="hover", options=list(container="body")),
         bsTooltip("tc_in_2", "The \"run-in\" (or \"wash-in\") is a measure of how long it takes for a treatment to ramp up to its full effect",
                   placement="right", trigger="hover", options=list(container="body"))
        ),
        conditionalPanel(
         condition = "input.n_treatments >= 3",
         hr(),
         h4("Treatment C"),
         numericInput("effect_3", "Effect size:", value=-40),
         numericInput("tc_in_3", HTML(paste0("Run-in ", tags$span("(days)"), ":")), min=1, value=3),
         numericInput("tc_out_3", HTML(paste0("Wash-out ", tags$span("(days)"), ":")), min=1, value=2),
         bsTooltip("effect_3", "The effect size is the asymptotic magnitude change from baseline that a person would experience in the long run, if they continued this treatment",
                   placement="right", trigger="hover", options=list(container="body")),
         bsTooltip("tc_out_3", "The \"wash-out\" (or \"carryover\") is a measure of how long it takes for a treatment to stop working once it is discontinued",
                   placement="right", trigger="hover", options=list(container="body")),
         bsTooltip("tc_in_3", "The \"run-in\" (or \"wash-in\") is a measure of how long it takes for a treatment to ramp up to its full effect",
                   placement="right", trigger="hover", options=list(container="body"))
        ),
        conditionalPanel(
         condition = "input.n_treatments >= 4",
         hr(),
         h4("Treatment D"),
         numericInput("effect_4", "Effect size:", value=-40),
         numericInput("tc_in_4", HTML(paste0("Run-in ", tags$span("(days)"), ":")), min=1, value=3),
         numericInput("tc_out_4", HTML(paste0("Wash-out ", tags$span("(days)"), ":")), min=1, value=2),
         bsTooltip("effect_4", "The effect size is the asymptotic magnitude change from baseline that a person would experience in the long run, if they continued this treatment",
                   placement="right", trigger="hover", options=list(container="body")),
         bsTooltip("tc_out_4", "The \"wash-out\" (or \"carryover\") is a measure of how long it takes for a treatment to stop working once it is discontinued",
                   placement="right", trigger="hover", options=list(container="body")),
         bsTooltip("tc_in_4", "The \"run-in\" (or \"wash-in\") is a measure of how long it takes for a treatment to ramp up to its full effect",
                   placement="right", trigger="hover", options=list(container="body"))
        ),
        conditionalPanel(
         condition = "input.n_treatments >= 5",
         hr(),
         h4("Treatment E"),
         numericInput("effect_5", "Effect size:", value=-40),
         numericInput("tc_in_5", HTML(paste0("Run-in ", tags$span("(days)"), ":")), min=1, value=3),
         numericInput("tc_out_5", HTML(paste0("Wash-out ", tags$span("(days)"), ":")), min=1, value=2),
         bsTooltip("effect_5", "The effect size is the asymptotic magnitude change from baseline that a person would experience in the long run, if they continued this treatment",
                   placement="right", trigger="hover", options=list(container="body")),
         bsTooltip("tc_out_5", "The \"wash-out\" (or \"carryover\") is a measure of how long it takes for a treatment to stop working once it is discontinued",
                   placement="right", trigger="hover", options=list(container="body")),
         bsTooltip("tc_in_5", "The \"run-in\" (or \"wash-in\") is a measure of how long it takes for a treatment to ramp up to its full effect",
                   placement="right", trigger="hover", options=list(container="body"))
        )
      ),
      column(4,
        h3("Study Design"),
        sliderInput("n_blocks", "Number of treatment blocks:", value=2, min=1, max=10),
        uiOutput("treatmentScheduleOptions"),
        numericInput("treatment_period", HTML(paste0("Treatment period ", tags$span("(days)"), ":")), value=30),
        numericInput("sampling_timestep", HTML(paste0("Sampling timestep ", tags$span("(days)"), ":")), min=0.1, value=1, step=0.1),
        # Tooltip for treatment_period is defined in renderUI using tipify, due to timing/rendering issue
        bsTooltip("treatment_period", "Number of days for each treatment period",
                  placement="right", trigger="hover", options=list(container="body")),
        bsTooltip("sampling_timestep", "Number of days between samples",
                  placement="right", trigger="hover", options=list(container="body")),
        hr(),
        h3("Noise Parameters"),
        numericInput("sd_baseline", HTML(paste0("Baseline drift ", tags$span("(SD)"), ":")), min=0, value=0.4, step=0.1),
        numericInput("sd_outcome", HTML(paste0("Process noise ", tags$span("(SD)"), ":")), min=0, value=0.6, step=0.1),
        numericInput("sd_obs", HTML(paste0("Observation noise ", tags$span("(SD)"), ":")), min=0, value=4, step = 0.1),
        bsTooltip("sd_baseline", "Baseline drift is a steady increase or decrease from baseline caused by long-term processes (e.g. long-term illness onset and recovery)",
                  placement="right", trigger="hover", options=list(container="body")),
        bsTooltip("sd_outcome", "Process noise consists of short-term fluctuations (e.g. changes in sleep and diet from day to day)",
                  placement="right", trigger="hover", options=list(container="body")),
        bsTooltip("sd_obs", "Observation noise is a function of the measurement instrument and not the underlying biological process (e.g. imprecision in a blood pressure cuff)",
                  placement="right", trigger="hover", options=list(container="body"))
      ),
      column(4,
        h3("Randomization"),
        numericInput("random_seed", "Random seed:", 1),
        actionButton("random_seed_generate", "Generate"),
        bsTooltip("random_seed", "The random seed used to generate the simulation. The choice of random seed can be used to exactly reproduce a simulation.",
                  placement="right", trigger="hover", options=list(container="body")),
        bsTooltip("random_seed_generate", "Generate a random seed",
                  placement="right", trigger="hover", options=list(container="body")),
        hr(),
        h3("Save data"),
        bookmarkButton(id="bookmark_button", label="Bookmark link", title="Bookmark this application's state and get a URL for sharing"),
        downloadButton("downloadData"),
        bsTooltip("bookmark_button", "",
                  placement="right", trigger="hover", options=list(container="body")),
        bsTooltip("downloadData", "Download timeseries data as a CSV file",
                  placement="right", trigger="hover", options=list(container="body")),
        hr(),
        actionButton("reset_button", "Reset inputs to default")
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
    input$sampling_timestep
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
      baseline_initial = 0,
      effect_size = effect_sizes(),
      tc_in = tc_ins(),
      tc_out = tc_outs(),
      tc_outcome = 1,
      sd_baseline = input$sd_baseline,
      sd_outcome = input$sd_outcome,
      sd_obs = input$sd_obs,
      treatment_period = input$treatment_period,
      sampling_timestep = input$sampling_timestep,
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
    tipify(selectInput("treatment_schedule", "Treatment schedule:", treatment_options()),
           title="Ordering of treatment periods", placement="right", trigger="hover", options=list(container="body"))
  })
  
  output$outcomePlot <- renderPlot({
    # prevent plotting before treatment_schedule has been set
    req(input$treatment_schedule)
    
    timeseries <- n1sim_result()$timeseries
    plot <- ggplot(data = timeseries, aes(x = t)) 
    plot <- plot + geom_line(aes(y = outcome, alpha="underlying"), colour="#555555", size=0.75)
    tmtfactor <- LETTERS[factor(timeseries$treatment)]
    plot <- plot + geom_point(aes(y = outcome_obs, alpha="observed", colour=tmtfactor, shape=tmtfactor, fill=tmtfactor))
    
    isolate({
      grid_x_major = seq(0, input$treatment_period*input$n_blocks*input$n_treatments, input$treatment_period)
      grid_x_minor = seq(input$treatment_period/2, input$treatment_period*(input$n_blocks*input$n_treatments-1/2), input$treatment_period)
    })
    plot <- plot + scale_x_continuous(minor_breaks = grid_x_minor, breaks = grid_x_major)
    plot <- plot + ggtitle("Underlying and observed outcome") + theme(plot.title = element_text(hjust=0.5, face="bold"))
    plot <- plot + xlab("Day in study") + ylab("Outcome")
    plot <- plot + theme(text = element_text(size=14))
    plot <- plot + scale_colour_manual(name="Treatment",
                                       values=c("A"="#009E73","B"="#D55E00","C"="#0072B2","D"="#CC79A7","E"="#E69F00"))
    plot <- plot + scale_fill_manual(name="Treatment",
                                     values=c("A"="#009E73","B"="#D55E00","C"="#0072B2","D"="#CC79A7","E"="#E69F00"))
    plot <- plot + scale_shape_manual(name="Treatment",
                                      values=c("A"=21,"B"=22,"C"=23,"D"=24,"E"=25))
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
    updateNumericInput(session, "sampling_timestep", value = 1)
    
    updateNumericInput(session, "sd_baseline", value = 0.4)
    updateNumericInput(session, "sd_outcome", value = 0.6)
    updateNumericInput(session, "sd_obs", value = 4)
    
    updateNumericInput(session, "random_seed", value = 1)
  })
  
  ## Bookmark state
  setBookmarkExclude("bookmark_button")
  observeEvent(input$bookmark_button, {
    session$doBookmark()
  })
  
  ## Downloadable csv of dataset 
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("n1sim-timeseries_random-seed-", input$random_seed, ".csv")
    },
    content = function(file) {
      result = n1sim_result()
      write.csv(result$timeseries, file)
      #save(result, file="n1sim-result.RData")
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server, enableBookmarking="url")

