rij = function(i, j, pi.clust, theta, theta.sd, theta.clust, junk.mixture, df, mu, sig, null.mixture, mu.null, sig.null){
  tmp.lgt = length(i);
  if(junk.mixture | null.mixture){
    tmp.lgt2 = length(j) - sum(junk.mixture) - sum(null.mixture);
    tmp.mn = theta.clust[1:tmp.lgt2];
    j = 1:(length(j) - sum(junk.mixture) - sum(null.mixture));
  }else{
    tmp.lgt2 = length(j);
    tmp.mn = theta.clust;
  }
  if(junk.mixture | null.mixture){
    tmp.lg.num = log.norm(i, j, tmp.lgt, tmp.lgt2, theta, theta.sd, tmp.mn);
    if(null.mixture){
      tmp.lg.num = c(tmp.lg.num, null.den(x = theta, mu = mu.null, sig = sig.null, log = TRUE));
    }
    if(junk.mixture){
      tmp.lg.num = c(tmp.lg.num, gen.t(theta, df = df, mu = mu, sig = sig, log = TRUE));
    }
    num = exp(matrix(rep(log(pi.clust), each = tmp.lgt) + tmp.lg.num, nrow = tmp.lgt, ncol=  tmp.lgt2 + sum(junk.mixture) + sum(null.mixture)));
    num[num <1e-300] = 1e-300;
    den = rowSums(num);
  }else{
    tmp.lg.num = log.norm(i, j = 1:length(theta.clust), tmp.lgt, tmp.lgt2 = length(theta.clust), theta, theta.sd, theta.clust);
    num = exp(matrix(rep(log(pi.clust[j]), each = tmp.lgt) + tmp.lg.num, nrow = tmp.lgt, ncol=  tmp.lgt2));
    num[num <1e-300] = 1e-300;
    den = rowSums(num);
  }
  res = num/den;
  return(res)
}
