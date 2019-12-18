pi.updt = function(j, m, pi.clust, theta, theta.sd, theta.clust, junk.mixture, df, mu, sig, junk.prob, fix.junk.prob,
                   null.mixture, mu.null, sig.null, null.prob, fix.null.prob){
  tmp = rij(1:m, j, pi.clust, theta, theta.sd, theta.clust, junk.mixture, df, mu, sig, null.mixture, mu.null, sig.null);
  num = if(length(j)>1){colSums(tmp)}else{sum(tmp)};
  jk.null = sum(junk.mixture) + sum(null.mixture);
  n.prob = j.prob = FALSE;
  if(null.mixture & !is.null(null.prob)){
    tmp.lgt = length(num) - sum(junk.mixture);
    if(num[tmp.lgt]/m < null.prob | fix.null.prob == TRUE){
      num[tmp.lgt] = null.prob*m;
      n.prob = TRUE;
    }
  }
  if(junk.mixture & !is.null(junk.prob)){
    tmp.lgt = length(num);
    if(num[length(num)]/m < junk.prob | fix.junk.prob == TRUE){
      num[tmp.lgt] = junk.prob*m;
      j.prob = TRUE;
    }
  }
  if(j.prob | n.prob){
    num[1:(length(num)-jk.null)] = m*(1 - sum(num[(length(num)- jk.null + 1):length(num)])/m)*num[1:(length(num)-jk.null)]/(sum(num[1:(length(num)-jk.null)]));
  }
  return(num/m);
}
