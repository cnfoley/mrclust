# introduce "scale" variable to shrink denisty to fit in the vertical plot range
scaled.dnorm <- function(theta, mean = 4, sd = 1, log = TRUE, scale = 1) {
  tmp <- stats::dnorm(theta, mean, sd)
  if (log) {
    tmp <- log(tmp)
  }
  return(scale * tmp)
}
