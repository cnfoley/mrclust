best.clust = function(dta){
  res = NULL;
  unique.obs = unique(dta$observation);
  for(i in unique.obs){
    tmp = dta$observation==i;
    mx = which.max(dta$probability[tmp]);
    res = rbind(res, dta[tmp,][mx,])
  }
  return(res);
}
