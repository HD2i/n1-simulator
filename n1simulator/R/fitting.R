library(nlme)

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
