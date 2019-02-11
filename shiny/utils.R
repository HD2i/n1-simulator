if(!require(arrangements)){
  install.packages("arrangements")
  library(arrangements)
}

generate_treatment_schedule_options <- function(n_treatments, n_blocks) {
  treatment_perms <- arrangements::permutations(n_treatments)
  order_perms <- arrangements::permutations(nrow(treatment_perms),n_blocks,replace=TRUE)
  order_perms_vec <- as.vector(t(order_perms))
  schedule_perms <- vector()
  for (i in order_perms_vec) {
    schedule_perms <- c(schedule_perms, treatment_perms[i,])
  }
  schedule_perms = matrix(schedule_perms,nrow=nrow(order_perms),byrow=TRUE)
}

treatment_schedule_options_to_strvec <- function(matrix) {
  options = vector()
  for (i in 1:nrow(matrix)) {
    entry <- setNames(paste(matrix[i,], collapse=""),paste(LETTERS[matrix[i,]], collapse=""))
    options <- c(options,entry)
  }
  return(options)
}