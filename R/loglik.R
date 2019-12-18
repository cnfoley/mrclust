loglik = function(m, pi.clust, theta, theta.sd, theta.clust, junk.mixture, df, mu, sig, null.mixture, mu.null, sig.null){
  i = 1:m;
  if(junk.mixture | null.mixture){
    K = length(theta.clust) - sum(junk.mixture) - sum(null.mixture);
  }else{
    K = length(theta.clust);
  }
  tmp.lg.den = log.norm(i, j = 1:K, m, tmp.lgt2 = K, theta, theta.sd, theta.clust);
  if(null.mixture){
    tmp.lg.den = c(tmp.lg.den, null.den(x = theta, mu = mu.null, sig = sig.null, log = TRUE));
  }
  if(junk.mixture){
    tmp.lg.den = c(tmp.lg.den, gen.t(theta, df = df, mu = 0, sig = sig, log = TRUE));
  }

  tmp = exp(matrix(rep(log(pi.clust), each  = m) + tmp.lg.den, nrow = m, ncol = length(pi.clust)));
  tmp.rowsum = rowSums(tmp)
  res = sum(log(tmp.rowsum));
  return(res)
}
