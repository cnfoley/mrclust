# introduce "scale" variable to shrink denisty to fit in the vertical plot range
gen_t_scale <- function(theta, df = 4, mu, sig, log = TRUE, scale = 1) {
  tmp <- ((gamma((df + 1) / 2) / gamma(df / 2) / sqrt(pi * df) / sig) *
            (1 + (theta - mu)^2 / sig^2 / df) ^ (- (df + 1) / 2))
  if (log) {
    tmp <- log(tmp)
  }
  return(scale * tmp)
}
