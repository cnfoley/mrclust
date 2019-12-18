null.den = function(x, mu, sig, log = TRUE){
  tmp = -0.5*log(2*pi) - log(sig) - 0.5*(x - mu)^2/sig^2;
  return(tmp)
}
