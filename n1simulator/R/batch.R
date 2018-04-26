n1_expand_parameters_and_run_experiment <- function(
  n_treatments,
  n_blocks,
  baseline_initial, effect_size,
  tc_in, tc_out, tc_outcome,
  sd_baseline, sd_outcome, sd_obs,
  treatment_period, sampling_timestep, noise_timestep,
  n_replicates,
  treatment_mat_by_block = NULL,
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
    n_replicates
  )
  n1_run_experiment(n_treatments, params, treatment_mat_by_block, initial_random_seed, baseline_func, cores)
}

n1_simulate_and_fit <- function(
  n_blocks, n_treatments,
  baseline_initial, effect_size_vec,
  tc_in_vec, tc_out_vec, tc_outcome,
  sd_baseline, sd_outcome, sd_obs,
  treatment_period, sampling_timestep, noise_timestep,
  treatment_mat_by_block = NULL,
  random_seed = NA,
  baseline_func = NULL
) {
  data <- n1_simulate(
    n_blocks, n_treatments, baseline_initial, effect_size_vec, tc_in_vec, tc_out_vec, tc_outcome,
    sd_baseline, sd_outcome, sd_obs,
    treatment_period, sampling_timestep, noise_timestep,
    treatment_mat_by_block,
    random_seed, baseline_func,
    return_data_frame = TRUE
  )

  fits <- lapply(1:4, function(i) n1_fit(n_treatments, data, i))
  list(
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

n1_expand_parameters <- function(
  n_treatments,
  n_blocks,
  baseline_initial, effect_size,
  tc_in, tc_out, tc_outcome,
  sd_baseline, sd_outcome, sd_obs,
  treatment_period, sampling_timestep, noise_timestep,
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
        replicate_id = 1:n_replicates
      )
    ))
  ))
  params
}

get_vec_columns <- function(df, prefix, n_cols) {
  as.matrix(df[sapply(1:n_cols, function(i) sprintf('%s_%d', prefix, i))])
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

n1_run_experiment <- function(n_treatments, parameters, treatment_mat_by_block = NULL, initial_random_seed = NA, baseline_func = NULL, cores = 1) {
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
  
  # TODO: support sweeping over treatment ordering
  treatment_order_arr <- aperm(simplify2array(
    lapply(1:n_trials, function(i) {
      if(is.null(treatment_mat_by_block)) {
        randomize_treatments_by_block(p$n_blocks[i], n_treatments)
      }
      else {
        treatment_mat_by_block
      }
    })
  ), c(3, 1, 2))
  
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
        treatment_mat_by_block = treatment_order_arr[i,,],
        random_seed = random_seeds[i],
        baseline_func = baseline_func
      )
    }
  )

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
      replicate_id = array(p$replicate_id),
      random_seed = array(random_seeds),
      treatment_order = treatment_order_arr
    ),
    combine_arrays(fit_results)
  )
}

