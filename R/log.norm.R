log.norm = function(i, j, tmp.lgt, tmp.lgt2, theta, theta.sd, theta.clust){
  tmp.theta = rep(theta[i], tmp.lgt2);
  tmp.theta.sd = rep(theta.sd[i], tmp.lgt2);
  tmp.clust = rep(theta.clust[j], each = tmp.lgt)
  res = -0.5*log(2*pi) - log(tmp.theta.sd) - 0.5*(tmp.theta - tmp.clust)^2/tmp.theta.sd^2;
  #mx = which.max(res);
  #return(list(exp(res - res[mx]), res[mx]));
  return(res);
}
