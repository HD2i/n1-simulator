# Script for generating simulated data for cognition study
# Based on https://github.com/HD2i/n1-simulator/examples/cognition-study/simulate-cognitive-study-data.Rmd by Noah Zimmerman
#
# Author: Mike Jones
# Date: 8-Sept-2019

library(ggplot2)
library(dplyr)
source('n1-simulator.R')

### Constants
N_BLOCKS = 2 # baseline included in block 1
N_TREATMENTS = 3 # includes baseline
TC_IN = c(0.1, 0.1, 0.1)
TC_OUT = c(0.1, 0.1, 0.1)
TC_OUTCOME = 0.1
SD_BASELINE = 0.05
#SD_OUTCOME = 0.2 # sd_outcome defined for each test below
SD_OBS = 0.1
SAMPLING_TIMESTEP = 1
NOISE_TIMESTEP = 0.1
TREATMENT_ORDER = rbind(c(1,2,3),c(1,3,2))
#RANDOM_SEED = 5 # defaults to a different random seed for each run
RETURN_DATA_FRAME = TRUE


### Post-processing utility functions
censor_treatment_period <- function(n1_simulate_result, treatment_period ){
  rows.to.censor = seq( from = 1, to = nrow(n1_simulate_result), 
                        by = treatment_period )
  slice(n1_simulate_result, -rows.to.censor)
}

censor_extra_points <- function(n1_simulate_result, initial_treatment_period, target_baseline_period, target_treatment_period){
  rows.to.keep = seq(from=2, to=1+target_baseline_period, by=1) # baseline treatment period
  n_treatment_periods = N_BLOCKS*N_TREATMENTS
  for(i in 1:(n_treatment_periods-1)) {
    rows.to.keep = c(rows.to.keep, seq(from=2+initial_treatment_period*i, to=1+initial_treatment_period*i+target_treatment_period, by=1))
  }
  slice(n1_simulate_result, rows.to.keep)
}

convert_treatment1_to_baseline <- function( n1_simulate_result ){
  block1 = filter(n1_simulate_result, block == 1)
  rest   = filter(n1_simulate_result, block != 1) %>%
    filter( treatment != 1)
  
  rbind(block1, rest) %>%
    mutate(time = row_number()) %>%
    select(time, block, treatment, outcome_obs, outcome)
}

lower_bound_data <- function(n1_simulate_result, lower_bound) {
  n1_simulate_result[['outcome']] <- replace(n1_simulate_result[['outcome']], n1_simulate_result[['outcome']] < lower_bound, lower_bound)
  n1_simulate_result[['outcome_obs']] <- replace(n1_simulate_result[['outcome_obs']], n1_simulate_result[['outcome_obs']] < lower_bound, lower_bound)
  return(n1_simulate_result)
}


### Data generator (wrapper for n1_simulate function with data post-processing)
generator <- function(initial_treatment_period, target_baseline_period, target_treatment_period, 
                      baseline_initial, effect_size, sd_outcome, lower_bound=NULL){
  
  simulated.data <- n1_simulate(
    n_treatments = N_TREATMENTS,
    n_blocks = N_BLOCKS,
    baseline_initial = baseline_initial,
    effect_size = effect_size, 
    tc_in = TC_IN,
    tc_out = TC_OUT,
    tc_outcome = TC_OUTCOME,
    sd_baseline = SD_BASELINE,
    sd_outcome = sd_outcome,
    sd_obs = SD_OBS,
    treatment_period = initial_treatment_period, 
    sampling_timestep = SAMPLING_TIMESTEP, 
    noise_timestep = NOISE_TIMESTEP,
    treatment_order = TREATMENT_ORDER, 
    return_data_frame = RETURN_DATA_FRAME
  )
  
  processed.data = censor_extra_points(simulated.data$timeseries,
                                       initial_treatment_period,
                                       target_baseline_period,
                                       target_treatment_period)
  processed.data = convert_treatment1_to_baseline(processed.data)
  
  if(!is.null(lower_bound)) {
    processed.data = lower_bound_data(processed.data, lower_bound)
  }
  
  processed.data
}
 
 
### Script parameters

# initial_treatment period defines number of days per treatment period to start with 
# Post-processing removes extra samples from each treatment period based on target treatment period lengths
# This should be greater than any number in target_baseline_period_vec or target_treatment_period_vec
initial_treatment_period = 10 

# 3 study durations (5, 15, or 27 days) and the corresponding length of baseline and treatment periods
study_duration_vec = c(5,15,27)
target_baseline_period_vec = c(1,3,7)
target_treatment_period_vec = c(1,3,5)

# This list parametrizes the script to generate data for different cognition tests
# Each test is comprised of multiple result types
# Each result type defines parameters like baseline_initial, and defines a list of scenarios for different effect size vectors
# Under each scenario, a description is included for quick reference. "A" = treatment A, "B" = treatment B, "N" = baseline.
test_definition_list = list(
  rat = list(
    score_percentage = list(
      baseline_initial = 35,
      sd_outcome = 7,
      lower_bound = 0,
      scenarios = list(
        A = list(
          effect_size = c(0, 55, 35),
          description = "A > N, B > N, A > B"
        ),
        B = list(
          effect_size = c(0, 35, -25),
          description = "A > N, B < N, A > B"
        ),
        C = list(
          effect_size = c(0, 35, 0),
          description = "A > N, B = N, A > B"
        ),
        D = list(
          effect_size = c(0, 35, 35),
          description = "A > N, B > N, A = B"
        ),
        E = list(
          effect_size = c(0, 0, 0),
          description = "A = N, B = N, A = B"
        )
      )
    ),
    avrg_response_time = list(
      baseline_initial = 18,
      sd_outcome = 2,
      lower_bound = 0,
      scenarios = list(
        A = list(
          effect_size = c(0, -10, -8),
          description = "A > N, B > N, A > B"
        ),
        B = list(
          effect_size = c(0, -8, 8),
          description = "A > N, B < N, A > B"
        ),
        C = list(
          effect_size = c(0, -8, 0),
          description = "A > N, B = N, A > B"
        ),
        D = list(
          effect_size = c(0, -8, -8),
          description = "A > N, B > N, A = B"
        ),
        E = list(
          effect_size = c(0, 0, 0),
          description = "A = N, B = N, A = B"
        )
      )
    )
  ),
  stroop = list(
    num_correct = list(
      baseline_initial = 5,
      sd_outcome = 0.5,
      lower_bound = 0,
      scenarios = list(
        A = list(
          effect_size = c(0, 4, 3),
          description = "A > N, B > N, A > B"
        ),
        B = list(
          effect_size = c(0, 3, -3),
          description = "A > N, B < N, A > B"
        ),
        C = list(
          effect_size = c(0, 3, 0),
          description = "A > N, B = N, A > B"
        ),
        D = list(
          effect_size = c(0, 3, 3),
          description = "A > N, B > N, A = B"
        ),
        E = list(
          effect_size = c(0, 0, 0),
          description = "A = N, B = N, A = B"
        )
      )
    ),
    total_time = list(
      baseline_initial = 18,
      sd_outcome = 1,
      lower_bound = 0,
      scenarios = list(
        A = list(
          effect_size = c(0, -12, -9),
          description = "A > N, B > N, A > B"
        ),
        B = list(
          effect_size = c(0, -9, 9),
          description = "A > N, B < N, A > B"
        ),
        C = list(
          effect_size = c(0, -9, 0),
          description = "A > N, B = N, A > B"
        ),
        D = list(
          effect_size = c(0, -9, -9),
          description = "A > N, B > N, A = B"
        ),
        E = list(
          effect_size = c(0, 0, 0),
          description = "A = N, B = N, A = B"
        )
      )
    )
  ),
  trailmaking = list(
    num_errors = list(
      baseline_initial = 5,
      sd_outcome = 1,
      lower_bound = 0,
      scenarios = list(
        A = list(
          effect_size = c(0, -4, -3),
          description = "A > N, B > N, A > B"
        ),
        B = list(
          effect_size = c(0, -3, 3),
          description = "A > N, B < N, A > B"
        ),
        C = list(
          effect_size = c(0, -3, 0),
          description = "A > N, B = N, A > B"
        ),
        D = list(
          effect_size = c(0, -3, -3),
          description = "A > N, B > N, A = B"
        ),
        E = list(
          effect_size = c(0, 0, 0),
          description = "A = N, B = N, A = B"
        )
      )
    ),
    total_time = list(
      baseline_initial = 15,
      sd_outcome = 1,
      lower_bound = 0,
      scenarios = list(
        A = list(
          effect_size = c(0, -10, -8),
          description = "A > N, B > N, A > B"
        ),
        B = list(
          effect_size = c(0, -8, 8),
          description = "A > N, B < N, A > B"
        ),
        C = list(
          effect_size = c(0, -8, 0),
          description = "A > N, B = N, A > B"
        ),
        D = list(
          effect_size = c(0, -8, -8),
          description = "A > N, B > N, A = B"
        ),
        E = list(
          effect_size = c(0, 0, 0),
          description = "A = N, B = N, A = B"
        )
      )
    )
  )
)

### Script
output_path = file.path(getwd(),"output")
dir.create(output_path)

for (i in 1:length(study_duration_vec)) {
  study_duration = study_duration_vec[i]
  target_baseline_period = target_baseline_period_vec[i]
  target_treatment_period = target_treatment_period_vec[i]
  
  for (test_name_i in names(test_definition_list)) {
    result_types = test_definition_list[[test_name_i]]
    
    for (result_type_i in names(result_types)) {
      result_type = result_types[[result_type_i]]
      
      for (scenario_i in names(result_type$scenarios)) {
        scenario = result_type$scenarios[[scenario_i]]
        
        output = generator(initial_treatment_period = initial_treatment_period,
                              target_baseline_period = target_baseline_period,
                              target_treatment_period = target_treatment_period,
                              baseline_initial = result_type$baseline_initial,
                              effect_size = scenario$effect_size,
                              sd_outcome = result_type$sd_outcome,
                              lower_bound = result_type$lower_bound
        )
        
        plot_title = paste(study_duration, "day", test_name_i, result_type_i, scenario$description)
        plot = ggplot(data = output, aes(x = time)) + 
          geom_line(aes(y = outcome)) + 
          geom_point(aes(y = outcome_obs, col = factor(treatment))) +
          ggtitle(plot_title) + 
          theme_bw() + theme(legend.position = "none")
        #print(plot)
        
        output_filename = paste(study_duration, "day", sep="")
        output_filename = paste(output_filename, test_name_i, result_type_i, scenario_i, sep="-")
        
        ggsave(paste0("plot-",output_filename,".png"), plot=plot, path=output_path)
        write.csv(output, file.path(output_path, paste0("output-",output_filename,".csv")))
      }
    }
  }
}

