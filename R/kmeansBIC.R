kmeansBIC = function(fit){

  m = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*log(m)*k)
}
