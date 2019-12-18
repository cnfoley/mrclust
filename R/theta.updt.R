theta.updt = function(j, m, theta, theta.sd, theta.clust, pi.clust, junk.mixture, df, mu, sig,
                      null.mixture, mu.null, sig.null){
  tmp = rij(1:m, j, pi.clust, theta, theta.sd, theta.clust, junk.mixture, df, mu, sig, null.mixture, mu.null, sig.null)/theta.sd[1:m]^2;
  den = if(length(j)>1){colSums(tmp)}else{sum(tmp)};
  num = if(length(j)>1){colSums(tmp*theta)}else{sum(tmp*theta)};
  null.jk = sum(null.mixture) + sum(junk.mixture);
  if(null.mixture & null.jk ==2){
    num[length(j) - 1] = mu.null;
  }else if(null.mixture & null.jk == 1){
    num[length(j)] = mu.null;
  }
  if(junk.mixture){num[length(j)] = mu}
  return(num/den);
}
