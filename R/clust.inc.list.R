clust.inc.list = function(dta, by.prob = 0.2, bound = 0){
  if(!((1-bound)/by.prob) %% 1==0){stop("(1-bound)/by.prob needs to be an integer");}
  unique.clust = unique(dta$cluster);
  res = vector('list', length(unique.clust));
  clust.nam = vector('character', length(unique.clust));
  nms.tmp = unique(dta$cluster.class);
  for(i in unique.clust){
    if(bound==0){sub.grps = (1/by.prob)}else{sub.grps = ceiling((1-bound)/by.prob)};
    tmp = vector('list', sub.grps);
    tmp.nam = vector('character', sub.grps);
    for(j in 1:sub.grps){
      p.mx = 1 - (j-1)*by.prob;
      p.mn = 1 - j*by.prob;
      if(p.mx<=bound){break}
      if((!bound==0) & (1 - j*by.prob < bound)){p.mn = bound;}
      tmp[[j]] = dta$observation[(dta$probability<=p.mx) & (dta$probability>p.mn) & (dta$cluster==i)];
      tmp.nam[j] = paste0("pr_range_",p.mx,"_to_",p.mn);
    }
    names(tmp) = tmp.nam;
    clust.nam[i] =  paste0("cluster_",nms.tmp[i]);
    res[[i]] = tmp;
  }
  names(res) = clust.nam;
  return(res)
}
