kmeansbic <- function(fit) {
  m <- length(fit$cluster)
  k <- nrow(fit$centers)
  d <- fit$tot_withinss
  return(d + 2 * log(m) * k)
}
