### SIMULATION ###

exponential_decay <- function(X_initial, X_target, tau, delta_t) {
  X_target + (X_initial - X_target) * exp(-delta_t / tau)
}

# Returns a vectorized function of time for the (deterministic) effect process
# for a single treatment
make_effect_function <- function(X_initial, t_change_vec, treatment_vec, E, tau0, tau1) {
  t_initial = t_change_vec[1]
  t_final = t_change_vec[length(t_change_vec)]

  n_periods <- length(t_change_vec) - 1
  stopifnot(n_periods >= 1)

  t_period_start <- t_change_vec[1:n_periods]
  t_period_end <- t_change_vec[2:(n_periods + 1)]

  # Construct a lookup table with effect values at the start of each time period
  # by analytically integrating the ODE
  X_vec <- numeric(n_periods)
  X_vec[1] <- X_initial
  if(n_periods > 1) {
    for(i in 1:(n_periods - 1)) {
      X_target <- if(treatment_vec[i]) E else 0
      tau <- if(treatment_vec[i]) tau1 else tau0
      X_vec[i+1] <- exponential_decay(X_vec[i], X_target, tau, t_period_end[i] - t_period_start[i])
    }
  }

  function(t_vec) {
    stopifnot(all(t_vec >= t_initial))
    stopifnot(all(t_vec <= t_final))

    # Create a vector of indices into t_change_vec where
    # period_indices_halfopen[i] is j if t_change_vec[j] <= t_vec[i] <  t_change_vec[j+1]
    # period_indices_halfopen[i] will be NA if t_vec[i] == t_final
    t_greq_start <- outer(t_vec, t_period_start, '>=')
    t_less_end <- outer(t_vec, t_period_end, '<')
    period_indices_halfopen = apply(t_greq_start & t_less_end, 1, function(v) match(TRUE, v))

    # Construct a fully non-NA vector of indices into t_change_vec
    # by assigning cases where t_vec[i] == t_final to the index (n_periods)
    period_indices <- ifelse(t_vec == t_final, n_periods, period_indices_halfopen)

    # Analytically integrate the ODE for the vector of time values,
    # using matched vectors of initial conditions and parameters for matched time periods
    X_target <- ifelse(treatment_vec[period_indices], E, 0)
    tau <- ifelse(treatment_vec[period_indices], tau1, tau0)
    exponential_decay(X_vec[period_indices], X_target, tau, t_vec - t_change_vec[period_indices])
  }
}


make_treatment_function <- function(t_change_vec, treatment_vec) {
  # Put an extra treatment on the end so length(t_change_vec) == length(treatment_vec)
  treatment_vec <- c(treatment_vec, treatment_vec[length(treatment_vec)])

  approxfun(t_change_vec, treatment_vec, method = 'constant')
}

make_brownian_baseline_function <- function(B_initial, noise_sd, t_initial, t_final, dt) {
  t_vec <- seq(t_initial, t_final, dt)
  B_vec <- c(B_initial, B_initial + cumsum(rnorm(length(t_vec) - 1, mean = 0, sd = noise_sd * sqrt(dt))))
  approxfun(t_vec, B_vec, method = 'linear')
}

#' Simulate N-of-1 study and return results as continuous functions of time.
#'
#' @param t_change_vec Sorted vector of simulation *change* times.
#' @export
#' @md
n1_simulate_functions <- function(
  t_change_vec,
  # t_vec[1] is the start time; t_vec[length(t_vec)] is the end time;
  # intermediate times correspond to treatment changes.
  B_func,       # Function that returns the baseline B at a particular time.
  E_vec,        # Vector of effect sizes for treatments. length(E_vec) is the number of treatments.
  T_mat,        # Matrix of treatments for each time period. ncol(T_mat) == nrow(E_vec). nrow(T_mat) == length(t_vec) - 1.
  # T_vec[i, j] == 1 (0) if treatment j is active (inactive) from time t_vec[i] to t_vec[i + 1].
  tau0_vec,     # Vector of time constants for deactivation of treatments
  tau1_vec,     # Vector of time constants for activation of treatments
  tc_outcome,   # Timescale of outcome variable response to effects
  sigma_Z,      # S.D. of process noise for outcome variable
  dt,           # Simulation timestep
  X_initial_vec = NA, # Initial values of effects. length(X_initial_vec) == length(E_vec)
  Z_initial = NA      # Initial value of outcome.
) {
  stopifnot(all(is.finite(t_change_vec)))
  stopifnot(all(t_change_vec == sort(t_change_vec)))
  stopifnot(all(is.finite(B_func(t_change_vec))))

  n_periods <- length(t_change_vec) - 1
  t_initial = t_change_vec[1]
  t_final = t_change_vec[n_periods + 1]

  stopifnot(all(is.finite(E_vec)))

  n_treatments = length(E_vec)

  stopifnot(all((T_mat == 0) | (T_mat == 1)))
  stopifnot(length(dt) == 1 && dt > 0)

  if(is.na(X_initial_vec)) {
    X_initial_vec <- rep(0, n_treatments)
  }
  else {
    stopifnot(length(X_initial_vec) == n_treatments)
    stopifnot(all(is.finite(X_initial_vec)))
  }

  if(is.na(Z_initial)) {
    Z_initial <- B_func(t_initial)
  }
  else {
    stopifnot(length(Z_initial_vec) == 1)
    stopifnot(is.finite(Z_initial_vec))
  }

  # Extend last treatment row to final time, so now nrow(T_mat) == length(t_vec)
  T_mat <- cbind(T_mat, T_mat[,ncol(T_mat)])

  # Actual simulation timesteps
  t_vec <- seq(t_initial, t_final, dt)
  n_t <- length(t_vec)
  if(t_vec[length(t_vec)] < t_final) {
    t_vec <- c(t_vec, t_final)
  }

  # Baseline evaluated at all simulation timesteps
  B_vec <- B_func(t_vec)

  # Treatment functions (just to return for convenience)
  T_funcs <- lapply(1:n_treatments, function(j) {
    make_treatment_function(t_change_vec, T_mat[1:n_periods,j])
  })

  # Deterministic effect process functions
  X_funcs <- lapply(1:n_treatments, function(i) {
    make_effect_function(X_initial_vec[i], t_change_vec, T_mat[,i], E_vec[i], tau0_vec[i], tau1_vec[i])
  })

  # Effects evaluated at all simulation timesteps (one effect per column)
  X_mat <- do.call(cbind, lapply(X_funcs, function(X_func) X_func(t_vec)))

  # Draw process noise for simulation steps
  Z_noise <- rnorm(n_t - 1, sd = sigma_Z * sqrt(dt))

  # Simulate outcome
  Z_vec <- numeric(n_t)
  Z_vec[1] <- Z_initial
  for(i in 1:(n_t - 1)) {
    Z_vec[i+1] <- exponential_decay(Z_vec[i], B_vec[i] + sum(X_mat[i,]), tc_outcome, dt) + Z_noise[i]
  }
  Z_func <- approxfun(t_vec, Z_vec, method = 'linear')

  list(
    B_func = B_func,
    T_funcs = T_funcs,
    X_funcs = X_funcs,
    Z_func = Z_func
  )
}

n1_observe_outcome <- function(Z_func, t_obs_vec, sigma_obs) {
  Z_func(t_obs_vec) + rnorm(length(t_obs_vec), sd = sigma_obs)
}

randomize_treatments_by_block <- function(n_blocks, n_treatments) {
  do.call(rbind, lapply(1:n_blocks, function(i) sample(1:n_treatments)))
}

treatment_mat_by_block_to_binary <- function(mat_by_block) {
  do.call(
    rbind,
    lapply(1:nrow(mat_by_block), function(i) diag(ncol(mat_by_block))[mat_by_block[i,],])
  )
}

treatment_mat_binary_to_vec <- function(mat) {
  apply(mat, 1, function(row) which(as.logical(row)))
}

treatment_order_mat_to_vec <- function(treatment_order_mat) {
  as.numeric(t(treatment_order_mat))
}

treatment_order_vec_to_str <- function(treatment_order) {
  paste(as.character(treatment_order), collapse = "")
}

treatment_order_mat_to_str <- function(treatment_order_mat) {
  treatment_order_vec_to_str(
    treatment_order_mat_to_vec(treatment_order_mat)
  )
}

treatment_order_str_to_vec <- function(treatment_order_str) {
  if(is.character(treatment_order_str)) {
    as.numeric(unlist(strsplit(treatment_order_str, split = "")))
  }
  else {
    NULL
  }
}

treatment_order_vec_to_mat <- function(n_blocks, treatment_order) {
  as.matrix(treatment_order, nrow=n_blocks, byrow=TRUE)
}

matrix_to_column_list <- function(mat, prefix) {
  setNames(
    lapply(1:ncol(mat), function(i) mat[,i]),
    sapply(1:ncol(mat), function(i) sprintf("%s_%d", prefix, i))
  )
}

n1_simulate <- function(
  n_treatments,
  n_blocks,

  baseline_initial, effect_size,
  tc_in, tc_out, tc_outcome,

  sd_baseline, sd_outcome, sd_obs,

  treatment_period,
  sampling_timestep,
  noise_timestep,
  
  treatment_order = NULL,
  random_seed = NA,

  baseline_func = NULL,
  
  return_data_frame = TRUE
) {
  if(is.na(random_seed)) {
    random_seed <- sample(2^31 - 1, 1)
  }
  set.seed(random_seed)

  study_duration <- n_blocks * n_treatments * treatment_period
  t_change_vec <- seq(0, study_duration, treatment_period)
  n_periods <- length(t_change_vec) - 1

  # Construct function to return block number as a function of time
  block_vec <- c(unlist(lapply(1:n_blocks, function(i) rep(i, n_treatments))), n_blocks)
  block_func <- approxfun(t_change_vec, block_vec, method = 'constant')
  
  # If NULL, randomize treatment order
  if(is.null(treatment_order)) {
    treatment_mat_by_block <- randomize_treatments_by_block(n_blocks, n_treatments)
  }
  # If a matrix, use it directly
  else if(is.matrix(treatment_order)){
    stopifnot(nrow(treatment_order) == n_blocks && ncol(treatment_order) == n_treatments)
    treatment_mat_by_block <- treatment_order
  }
  else {
    stopifnot(is.numeric(treatment_order))
    # If a one-block vector, repeat it
    if(length(treatment_order) == n_treatments) {
      treatment_mat_by_block <- matrix(rep(treatment_order, n_blocks), nrow = n_blocks, byrow = TRUE)
    }
    # If an all-blocks vector, convert it to a matrix
    else {
      stopifnot(length(treatment_order) == n_treatments * n_blocks)
      treatment_mat_by_block <- matrix(treatment_order, nrow = n_blocks, byrow = TRUE)
    }
  }
  treatment_mat <- treatment_mat_by_block_to_binary(treatment_mat_by_block)

  if(is.null(baseline_func)) {
    baseline_func <- make_brownian_baseline_function(baseline_initial, sd_baseline, 0, study_duration, noise_timestep)
  }

  result <- n1_simulate_functions(
    t_change_vec, baseline_func, effect_size, treatment_mat,
    tc_in, tc_out, tc_outcome,
    sd_outcome, noise_timestep
  )

  t_obs_vec <- seq(0, study_duration, sampling_timestep)

  block_vec_t_obs <- block_func(t_obs_vec)

  baseline_vec_t_obs <- baseline_func(t_obs_vec)

  treatment_mat_t_obs_binary <- do.call(cbind, lapply(result$T_funcs, function(T_func) T_func(t_obs_vec)))
  treatment_vec_t_obs <- treatment_mat_binary_to_vec(treatment_mat_t_obs_binary)

  effect_mat <- do.call(cbind, lapply(result$X_funcs, function(X_func) X_func(t_obs_vec)))
  effect_vec <- apply(effect_mat, 1, sum)

  outcome_obs_vec <- n1_observe_outcome(result$Z_func, t_obs_vec, sd_obs)

  outcome_vec <- result$Z_func(t_obs_vec)
  timeseries_obj <- if(return_data_frame) {
    # Construct data frame with a separate effect column for each treatment
    do.call(data.frame, c(
      list(
        t = t_obs_vec,
        block = block_vec_t_obs,
        baseline = baseline_vec_t_obs,
        treatment = treatment_vec_t_obs
      ),
      matrix_to_column_list(treatment_mat_t_obs_binary, "treatment"),
      list(
        effect = effect_vec
      ),
      matrix_to_column_list(effect_mat, "effect"),
      list(
        outcome = outcome_vec,
        outcome_obs = outcome_obs_vec
      )
    ))
  }
  else {
    # Construct a named list of time series, but with
    # the effect columns as a matrix `effect_by_treatment` with length(t_obs_vec) rows and n_treatments columns
    list(
      t = t_obs_vec,
      block = block_vec_t_obs,
      baseline = baseline_vec_t_obs,
      treatment = treatment_vec_t_obs,
      treatment_by_treatment = treatment_mat_t_obs_binary,
      effect = effect_vec,
      effect_by_treatment = effect_mat,
      outcome = outcome_vec,
      outcome_obs = outcome_obs_vec
    )
  }
  
  list(
    random_seed = random_seed,
    treatment_order = treatment_mat_by_block,
    timeseries = timeseries_obj
  )
}

treatment_mat_by_block_to_string <- function(treatment_mat_by_block, delimiter = "") {
  paste(as.character(t(treatment_mat_by_block)), collapse = delimiter)
}



### FITTING ###

n1_fit <- function(n_treatments, data, model) {
  if(model == 1) {
    n1_fit_model1(n_treatments, data)
  }
  else if(model == 2) {
    n1_fit_model2(n_treatments, data)
  }
  else if(model == 3) {
    n1_fit_model3(n_treatments, data)
  }
  else if(model == 4) {
    n1_fit_model4(n_treatments, data)
  }
  else if(is.function(model)) {
    # User-provided function
    model(n_treatments, data)
  }
  else {
    NULL
  }
}

n1_fit_model1 <- function(n_treatments, data) {
  format_lm_results(n_treatments, lm(outcome_obs ~ factor(treatment), data = data))
}

n1_fit_model2 <- function(n_treatments, data) {
  if(max(data$block) == 1) {
    format_lm_results(n_treatments, lm(outcome_obs ~ factor(treatment), data = data))
  }
  else {
    format_lm_results(n_treatments, lm(outcome_obs ~ factor(treatment) + factor(block), data = data))
  }
}

n1_fit_model3 <- function(n_treatments, data) {
  if(max(data$block) == 1) {
    format_lm_results(n_treatments, lm(outcome_obs ~ factor(treatment) + t * factor(treatment), data = data))
  }
  else {
    format_lm_results(n_treatments, lm(outcome_obs ~ factor(treatment) + t * factor(treatment) + factor(block), data = data))
  }
}

n1_fit_model4 <- function(n_treatments, data) {
  if(max(data$block) == 1) {
    format_lm_results(n_treatments, lm(outcome_obs ~ factor(treatment), data = data))
  }
  else {
    library(nlme)
    format_lme_results(n_treatments, lme(fixed = (outcome_obs ~ factor(treatment)), random = (~ 1 | block), data = data))
  }
}

format_lm_results <- function(n_treatments, fit_obj) {
  coeff_summ <- summary(fit_obj)$coefficients
  list(
    coefficients = data.frame(
      estimate = coeff_summ[1:n_treatments,1],
      pvalue = coeff_summ[1:n_treatments,4]
    ),
    fit_obj = fit_obj
  )
}

format_lme_results <- function(n_treatments, fit_obj) {
  coeff_summ <- summary(fit_obj)$tTable
  list(
    coefficients = data.frame(
      estimate = coeff_summ[1:n_treatments,1],
      pvalue = coeff_summ[1:n_treatments,5]
    ),
    fit_obj = fit_obj
  )
}



### BATCH MODE ###

n1_expand_parameters_and_run_experiment <- function(
  n_treatments,
  n_blocks,
  baseline_initial, effect_size,
  tc_in, tc_out, tc_outcome,
  sd_baseline, sd_outcome, sd_obs,
  treatment_period, sampling_timestep, noise_timestep,
  treatment_order = NULL,
  n_replicates,
  initial_random_seed = NA, baseline_func = NULL,
  cores = 1
) {
  params <- n1_expand_parameters(
    n_treatments,
    n_blocks,
    baseline_initial, effect_size,
    tc_in, tc_out, tc_outcome,
    sd_baseline, sd_outcome, sd_obs,
    treatment_period, sampling_timestep, noise_timestep,
    treatment_order,
    n_replicates
  )
  n1_run_experiment(n_treatments, params, initial_random_seed, baseline_func, cores)
}

n1_simulate_and_fit <- function(
  n_treatments, n_blocks,
  baseline_initial, effect_size,
  tc_in, tc_out, tc_outcome,
  sd_baseline, sd_outcome, sd_obs,
  treatment_period, sampling_timestep, noise_timestep,
  treatment_order = NULL,
  random_seed = NA,
  baseline_func = NULL
) {
  result <- n1_simulate(
    n_treatments, n_blocks, baseline_initial, effect_size, tc_in, tc_out, tc_outcome,
    sd_baseline, sd_outcome, sd_obs,
    treatment_period, sampling_timestep, noise_timestep,
    treatment_order,
    random_seed, baseline_func,
    return_data_frame = TRUE
  )
  data <- result$timeseries
  
  fits <- lapply(1:4, function(i) n1_fit(n_treatments, data, i))
  list(
    treatment_order = result$treatment_order,
    estimates = do.call(rbind, lapply(fits, function(fit) fit$coefficients[,'estimate'])),
    pvalues = do.call(rbind, lapply(fits, function(fit) fit$coefficients[,'pvalue']))
  )
}

paramlist_to_expandargs <- function(paramlist, prefix) {
  setNames(
    as.list(paramlist),
    sapply(1:length(paramlist), function(i) sprintf("%s_%d", prefix, i))
  )
}

expand_treatment_order <- function(treatment_order) {
  if(is.null(treatment_order)) {
    NA
  }
  else {
    if(is.list(treatment_order)) {
      unlist(lapply(treatment_order, function(toi) {
        if(is.vector(treatment_order)) {
          treatment_order_vec_to_str(toi)
        }
        else {
          treatment_order_mat_to_str(toi)
        }
      }))
    }
    else {
      if(is.vector(treatment_order)) {
        treatment_order_vec_to_str(treatment_order)
      }
      else {
        treatment_order_mat_to_str(treatment_order)
      }
    }
  }
}

n1_expand_parameters <- function(
  n_treatments,
  n_blocks,
  baseline_initial, effect_size,
  tc_in, tc_out, tc_outcome,
  sd_baseline, sd_outcome, sd_obs,
  treatment_period, sampling_timestep, noise_timestep,
  treatment_order = NULL,
  n_replicates
) {
  stopifnot(length(n_treatments) == 1)
  
  stopifnot(length(effect_size) == n_treatments)
  stopifnot(length(tc_in) == n_treatments)
  stopifnot(length(tc_out) == n_treatments)
  
  # Produce all combinations of parameters
  params <- rev(do.call(
    expand.grid,
    rev(c(
      stringsAsFactors = FALSE,
      list(
        n_blocks = n_blocks,
        baseline_initial = baseline_initial
      ),
      paramlist_to_expandargs(effect_size, 'effect_size'),
      paramlist_to_expandargs(tc_in, 'tc_in'),
      paramlist_to_expandargs(tc_out, 'tc_out'),
      list(
        tc_outcome = tc_outcome,
        sd_baseline = sd_baseline,
        sd_obs = sd_obs,
        sd_outcome = sd_outcome,
        treatment_period = treatment_period,
        sampling_timestep = sampling_timestep,
        noise_timestep = noise_timestep,
        treatment_order = expand_treatment_order(treatment_order),
        replicate_id = 1:n_replicates
      )
    ))
  ))
  params
}

get_vec_columns <- function(df, prefix, n_cols) {
  as.matrix(df[sapply(1:n_cols, function(i) sprintf('%s_%d', prefix, i))])
}

# Combines a list of arrays into a matrix with an extra dimension equal to the length fo the outer list.
listofarrays_to_array <- function(la) {
  stopifnot(is.list(la))
  
  # Get the initial and final dimensinos
  dim_initial <- dim(la[[1]])
  if(is.null(dim_initial)) {
    dim_combined <- length(la)
  }
  else {
    dim_combined <- c(dim_initial, length(la))
  }
  
  # Concatenate inner arrays using c()
  bigv <- do.call(c, la)
  
  # Turn it into an array with the right final dimensions
  biga <- array(bigv, dim = dim_combined)
  
  # Return the big array with the last dimension turned into the first
  aperm(biga, if(ndims == 1) 1 else c(ndims, 1:(ndims - 1)))
}

# Combines a list of list of arrays into a list of arrays,
# each of which has an extra initial dimension equal to the length of the original outer list.
#
# E.g.: list(list(n_blocks = 1, effect_size = array(c(2,3))), list(n_blocks = 2, effect_size = array(c(4, 5))))
# becomes something like list(n_blocks = array(c(1, 2)), effect_size = matrix(c(2, 3, 4, 5), byrow=TRUE))
#
# Length-1 numeric vectors are interpreted as 0-dimensional arrays (scalars).
# Longer numeric vectors are not allowed: they must be 1-dimensional arrays with explicit dimensions
# (to prevent ambiguity for length-1 numeric vectors).
combine_arrays  <- function(lla)  {
  stopifnot(is.list(lla))
  stopifnot(is.list(lla[[1]]))
  
  # Get the initial and final dimensinos
  dims_initial <- lapply(lla[[1]], function(a) dim(a))
  dims_combined <- lapply(dims_initial, function(dim_initial) {
    if(is.null(dim_initial)) {
      length(lla)
    }
    else {
      c(dim_initial, length(lla))
    }
  })
  
  # Concatenate inner arrays using c()
  l_bigv <- do.call(mapply, c(list(c), lla, list(SIMPLIFY = FALSE)))
  
  # Turn them into arrays with the right final dimensions
  l_biga <- mapply(
    function(bigv, dim_combined) {
      array(bigv, dim = dim_combined)
    },
    l_bigv, dims_combined,
    SIMPLIFY = FALSE
  )
  
  # Return list of arrays with original names and last dimension turned into the first
  lapply(l_biga, function(a) {
    ndims <- length(dim(a))
    aperm(a, if(ndims == 1) 1 else c(ndims, 1:(ndims - 1)))
  })
}

n1_run_experiment <- function(
  n_treatments, parameters,
  initial_random_seed = NA, baseline_func = NULL, cores = 1
) {
  if(is.na(initial_random_seed)) {
    initial_random_seed <- sample(2^31 - 1, 1)
  }
  
  p <- parameters
  n_trials <- nrow(p)
  trial_ids <- 1:n_trials
  random_seeds <- initial_random_seed:(initial_random_seed + n_trials - 1)
  
  effect_size_mat <- get_vec_columns(p, 'effect_size', n_treatments)
  tc_in_mat <- get_vec_columns(p, 'tc_in', n_treatments)
  tc_out_mat <- get_vec_columns(p, 'tc_out', n_treatments)
  
  if(cores > 1) {
    library(parallel)
    lapply_func <- mclapply
  }
  else {
    lapply_func <- lapply
  }
  
  fit_results <- lapply_func(
    trial_ids,
    function(i) {
      fit_result <- n1_simulate_and_fit(
        p$n_blocks[i], n_treatments,
        p$baseline_initial[i], effect_size_mat[i,],
        tc_in_mat[i,], tc_out_mat[i,], p$tc_outcome[i],
        p$sd_baseline[i], p$sd_outcome[i], p$sd_obs[i],
        p$treatment_period[i], p$sampling_timestep[i], p$noise_timestep[i],
        treatment_order = treatment_order_str_to_vec(p$treatment_order[[i]]),
        random_seed = random_seeds[i],
        baseline_func = baseline_func
      )
    }
  )
  #treatment_order <- array(unlist(lapply(fit_results, function(fr) fr$treatment_order)))
  #estimates <- listofarrays_to_array(lapply(fit_results, function(fr) fr$estimates))
  #pvalues <- listofarrays_to_array(lapply(fit_results, function(fr) fr$pvalues))
  
  c(
    list(
      trial_id = array(trial_ids),
      n_blocks = array(p$n_blocks),
      baseline_initial = array(p$baseline_initial),
      effect_size = effect_size_mat,
      tc_in = tc_in_mat,
      tc_out = tc_out_mat,
      tc_outcome = array(p$tc_outcome),
      sd_baseline = array(p$sd_baseline),
      sd_outcome = array(p$sd_outcome),
      sd_obs = array(p$sd_obs),
      treatment_period = array(p$treatment_period),
      sampling_timestep = array(p$sampling_timestep),
      noise_timestep = array(p$noise_timestep),
      #treatment_order = treatment_order,
      #estimate = estimates,
      #pvalue = pvalues,
      replicate_id = array(p$replicate_id),
      random_seed = array(random_seeds)
    ),
    combine_arrays(fit_results)
  )
}
