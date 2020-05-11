#' MR-Clust mixture model fitting
#'
#' Assessment of clustered heterogeneity in Mendelian randomization analyses using expectation-maximisation (EM) based model fitting of the MR-Clust mixture model. Function output includes both data-tables and a visualisation of the assingment of variants to clusters.
#'
#' @param theta numeric vector of length the number of variants, the i-th element is a ratio-estimate for the i-th genetic variant.
#' @param theta.se numeric vector of length the number of variants, the i-th element is the standard error of the ratio-estimate for the i-th genetic variant.
#' @param bx numeric vector of length the number of variants, the i-th element is the estimated regression coefficient - i.e. beta-x value - relating the i-th genetic variant to the risk-factor.
#' @param by numeric vector of length the number of variants, the i-th element is the estimated regression coefficient - i.e. beta-y value - relating the i-th genetic variant to the outcome.
#' @param bxse numeric vector of length the number of variants, the i-th element is the standard error of the estimated regression coefficient relating the i-th genetic variant to the risk-factor.
#' @param byse numeric vector of length the number of variants, the i-th element is the standard error of the estimated regression coefficient relating the i-th genetic variant to the outcome.
#' @param obs.names character vector of length the number of variants, the i-th element is the name of the i-th genetic variants - e.g. the rsID.
#' @param max.iter numeric integer denoting the maximum number of iterations to take before stopping the EM-algorithm's search for a maxima in the log-likelihood.
#' @param tol numeric scalar denoting the maximum absolute difference between two computations of the log-likelihood with which we accept that a maxima in the log-likelihood has been computed.
#' @param junk.sd numeric scalar denoting the scale parameter in the generalised t-distribution
#' @param junk.mean numeric scalar denoting the mean of the generalised t-distribution. By default mean is set to zero.
#' @param min.clust.search numeric integer which denotes the minimum number of clusters searched for in the data - default computes evidence supporting up to K=10 clusters which might explain any clustered heterogeneity in the data.
#' @param stop.bic.iter numeric integer I, for computational efficiency - particularly when analysing large numbers of variants - we can stop the EM-algorithm if the BIC is monotonic increasing over the previous I increases in the number of clusters K. By default evidence supporting at least 10 clusters in the data is computed and so, for example, if the BIC from models which assume 6 clusters; 7 clusters; ... or; 10 clusters is monotonic increasing - in the number of clusters K -then the EM-algorithm is stopped and the model whose K minimises the BIC is returned.
#' @param results.list character list allowing users to choose whether to return a table with the variants assigned to: "all" of the clusters; a single "best" cluster or; both. By default we return both, i.e. results.list = list("all", "best").
#' @param cluster.membership numeric list which allows users to output a list which, for each cluster, returns the variants assigned to the cluster by stratified by the probability of belonging to the cluster. By default, cluster.membership = list(by.prob = 0.1, bound = 0); so that MRClust returns a list, which for each cluster, outputs the variants assigned to the cluster with probability between (0.9,1); (0.8,0.9);... and finally; (0.1,0), i.e. by probability increments 0.1 from 1 to a lower bound of 0.
#' @param plot.results numeric list which allows users to plot the output of MRClust. By default, plot.results = list("best", min.pr = 0.5); so that the best clustering is plotted with variants assigned to a cluster with probability above 0.5.
#' @param trait.search logical, for each of the non-null and non-junk clusters search phenoscanner for traits associated with the variants.
#' @param trait.pvalue numeric scalar for use with trait.search, representing the maximum p-value with with at least one variant in the cluster must be associated with a trait for it to be returned in the phenoscanner search. Deafult value is GWA significance, i.e. 5*10^-8.
#' @param proxy.r2 numeric scalar for use with trait search, allowing variants whose r2>=proxy.r2 to be included in the trait search. Default r2=0.8.
#' @param catalogue character, for use with trait search. From Phenoscanner (http://www.phenoscanner.medschl.cam.ac.uk/information/) "the catalogue to be searched (options: None, GWAS, eQTL, pQTl, mQTL, methQTL)". Default setting is catalogue = "GWAS".
#' @param proxies character, for use with trait search. From Phenoscanner (http://www.phenoscanner.medschl.cam.ac.uk/information/) "the proxies database to be searched (options: None, AFR, AMR, EAS, EUR, SAS)". Default setting is proxies = "None"
#' @param build integer, for use with trait search. From Phenoscanner (http://www.phenoscanner.medschl.cam.ac.uk/information/) "Human genome build numbers (options: 37, 38; default: 37)". Default setting is build = 37.
#' @return Returned are: estimates of the putative number of clusters in the sample, complete with allocation probabilities and summaries of the association estimates for each variant; plots which visualise the allocation of variants to clusters and; several summaries of the fitting process, i.e. the BIC and likelihood estimates.
#' @export
mr_clust_em = function(theta, theta.se, bx, by, bxse, byse,
                       obs.names = NULL, max.iter = 5e3, tol = 1e-5,
                       junk.sd = NULL, junk.mean = 0,
                       stop.bic.iter = 5, min.clust.search = 10,
                       results.list = list("all", "best"),
                       cluster.membership = list(by.prob = 0.1, bound = 0),
                       plot.results = list("best", min.pr = 0.5),
                       trait.search = FALSE, trait.pvalue = 1e-5, proxy.r2 = 0.8,
                       catalogue = "GWAS", proxies = "None", build = 37){

  # previous user defined options now hard-coded
    cluster.sizes = NULL; init.clust.means = NULL; init.clust.probs = NULL;
    junk.mixture = TRUE; df = 4; junk.prob = NULL; junk.null.aware.bic = TRUE;
    fix.junk.prob = FALSE; null.mixture = TRUE; null.sd = NULL; null.prob = NULL; fix.null.prob = FALSE;
    scale.grid.search = FALSE; grid.increment = 0.1; grid.max = 2;
    clust.size.prior = FALSE; bic.prior = NULL;
    rand.init = TRUE; rand.num = 5; rand.sample = seq(0.05,0.4, by = 0.05);

  # initial parameter values

  m = length(theta);
  if(is.null(cluster.sizes)){cluster.sizes = 0:m};
  if(is.null(obs.names)){obs.names = paste0("snp_",1:m)};

  num.clust = length(cluster.sizes);
  if(rand.init){num.clust = num.clust*rand.num};

  bic.clust = bic.clust.mx = vector();
  clust.mn = vector('list', num.clust);
  clust.pi = vector('list', num.clust);
  log.like = vector('list', num.clust);


  count0 = count.K = 1;

  for(i in cluster.sizes){
    for(itr in 1:(rand.num+1)){

      if(is.null(init.clust.means) | is.null(init.clust.probs)){
        if(i > 0 & i != m){
          init.conds = kmeans(x = theta, centers = i, iter.max = 5e3);
          clust.means =  as.numeric(init.conds$centers);
          clust.probs = table(init.conds$cluster)/m;
          init.conds;clust.means;clust.probs
        }else if(i == m){
          clust.means =  theta;
          clust.probs = rep(1,m)/m;
        }else if(i ==0){
          clust.means = clust.probs = NULL;
        }
      }

      # junk distribution parameters

      if( junk.mixture & is.null(junk.sd)){
        rng.thet = range(theta)
        max.disp = which.max(abs(theta)+2*theta.se);
        sig = (rng.thet[2]-rng.thet[1] + theta.se[max.disp]);
      }else if( junk.mixture & !is.null(junk.sd)){
        sig = junk.sd;
      }else{
        sig = NULL;
      }
      if( junk.mixture & is.null(junk.mean)){
        mu = 0;
      }else if( junk.mixture & !is.null(junk.mean)){
        mu = junk.mean;
      }else{
        mu = NULL;
      }
      # null distribution parameters
      null.obs = abs(theta/theta.se) < 1.96;
      if(null.mixture & is.null(null.sd)){
        sig.null = theta.se;
        mu.null = 0;
      }else if( null.mixture & !is.null(null.sd)){
        sig.null = null.sd;
        mu.null = 0;
      }else{
        sig.null = mu.null =  NULL;
      }

      ####
      if(itr == 1){
        if(junk.mixture | null.mixture){
          ### junk.mixture parameters
          sum.pr = 0;
          junk.obs = sum(2*(1- pnorm(theta - median(theta), 0,theta.se))<0.05);
          #junk.obs = abs(theta/theta.se) < 1.96;
          if(junk.mixture & !fix.junk.prob){
            if(sum(junk.obs) == 0){
              jk.pr = 1/m;
              sum.pr = jk.pr;
            }else{
              jk.pr = sum(junk.obs)/m;
              sum.pr = jk.pr;
            }
          }else if(junk.mixture & fix.junk.prob){
            jk.pr = junk.prob;
            sum.pr = jk.pr;
          }else{
            jk.pr = NULL;
          }

          #### null.mixture parameters
          if(null.mixture & !fix.null.prob){
            if(is.null(jk.pr)){
              tmp.jk.pr = 0;
            }else{tmp.jk.pr = jk.pr;}
            if(sum(null.obs) == 0){
              null.pr = 1/m;
              sum.pr = sum.pr + null.pr;
            }else{
              null.pr = sum(null.obs)/m;
              if(null.pr + tmp.jk.pr > (1 - 1/m)){
                tmp.sum = null.pr + tmp.jk.pr;
                if(is.null(jk.pr)){
                  null.pr = null.pr/tmp.sum - 1/m;
                }else{
                  null.pr = null.pr/tmp.sum - 1/(2*m);
                  tmp.jk.pr = tmp.jk.pr/tmp.sum - 1/(2*m);
                  jk.pr = tmp.jk.pr;
                  sum.pr = tmp.jk.pr;
                }
              }
              sum.pr = sum.pr + null.pr;
            }
          }else if(null.mixture & fix.null.prob){
            null.pr = null.prob;
            sum.pr = sum.pr + null.pr;
          }else{
            null.pr = NULL;
          }

          clust.means =  c(clust.means, mu.null, mu);
          clust.probs =  c(clust.probs*(1 - sum.pr), null.pr, jk.pr);
          K.clust = cluster.sizes[count.K] + sum(junk.mixture) + sum(null.mixture);

        }else{
          K.clust = cluster.sizes[count.K];
          sig = mu = NULL;
          sig.null = mu.null = NULL;
          junk.null.aware.bic = FALSE;
        }
      }else{

        if(null.mixture){null.pr = sample(rand.sample, 1); sum.pr = null.pr;}else{null.pr = NULL; sum.pr = 0;}
        if(junk.mixture){jk.pr = sample(rand.sample, 1); sum.pr = sum.pr + jk.pr;}else{jk.pr = NULL;}

        clust.means =  c(clust.means, mu.null, mu);
        clust.probs =  c(clust.probs*(1 - sum.pr), null.pr, jk.pr);
        K.clust = cluster.sizes[count.K] + sum(junk.mixture) + sum(null.mixture);

      }
      ####

      if(scale.grid.search & junk.mixture){

        sigs = sig*seq(1, grid.max, by = grid.increment);
        pi.tmp = vector('list', length(sigs));
        theta.tmp = vector('list', length(sigs));
        loglik.tmp = vector('numeric', length(sigs));
        count2 = 1;

        for(inc in sigs){
          sig.tmp = inc;
          pi.clust = clust.probs;                # initialise cluster class probabilities
          theta.clust = clust.means;              # initialise cluster centroids; these can be the same as the std error is assumed to vary between the observations
          count = 1;                                        # start iteration count
          loglik.diff = 1;
          loglik.iter = NULL;
          while(loglik.diff > tol & count < max.iter){

            tmp.loglik0 = loglik(m, pi.clust, theta, theta.se, theta.clust, junk.mixture, df, mu, sig.tmp,
                                 null.mixture, mu.null, sig.null);
            tmp.clust = theta.updt(1:K.clust, m, theta, theta.se, theta.clust, pi.clust, junk.mixture, df, mu, sig.tmp,
                                   null.mixture, mu.null, sig.null);
            tmp.pi = pi.updt(1:K.clust, m, pi.clust, theta, theta.se, theta.clust, junk.mixture, df, mu, sig.tmp, junk.prob, fix.junk.prob,
                             null.mixture, mu.null, sig.null, null.prob, fix.null.prob);
            theta.clust = tmp.clust;
            pi.clust = tmp.pi;
            tmp.loglik = loglik(m, pi.clust, theta, theta.se, theta.clust, junk.mixture, df, mu, sig.tmp,
                                null.mixture, mu.null, sig.null);
            loglik.diff = abs(tmp.loglik - tmp.loglik0);
            loglik.iter = c(loglik.iter, tmp.loglik0);
            count = count + 1;
          }
          pi.tmp[[count2]] = pi.clust;
          theta.tmp[[count2]] = theta.clust;
          loglik.tmp[count2] = tmp.loglik;
          count2 = count2 + 1;
        }
        mx = which.max(loglik.tmp);
        theta.clust = theta.tmp[[mx]];
        pi.clust = pi.tmp[[mx]];
        tmp.loglik = loglik.tmp[mx];
        sig = sigs[mx];

      }else{

        pi.clust = clust.probs;                # initialise cluster class probabilities
        theta.clust = clust.means;              # initialise cluster centroids; these can be the same as the std error is assumed to vary between the observations
        count = 1;                                        # start iteration count
        loglik.diff = 1;
        loglik.iter = NULL;

        while(loglik.diff > tol & count < max.iter){

          tmp.loglik0 = loglik(m, pi.clust, theta, theta.se, theta.clust, junk.mixture, df, mu, sig, null.mixture, mu.null, sig.null);
          tmp.clust = theta.updt(1:K.clust, m, theta, theta.se, theta.clust, pi.clust, junk.mixture, df, mu, sig,
                                 null.mixture, mu.null, sig.null);
          tmp.pi = pi.updt(1:K.clust, m, pi.clust, theta, theta.se, theta.clust, junk.mixture, df, mu, sig, junk.prob, fix.junk.prob,
                           null.mixture, mu.null, sig.null, null.prob, fix.null.prob);
          theta.clust = tmp.clust;
          pi.clust = tmp.pi;
          tmp.loglik = loglik(m, pi.clust, theta, theta.se, theta.clust, junk.mixture, df, mu, sig, null.mixture, mu.null, sig.null);
          loglik.diff = abs(tmp.loglik - tmp.loglik0);
          loglik.iter = c(loglik.iter, tmp.loglik0);
          count = count + 1;
        }
      }

      clust.mn[[count0]] = theta.clust;
      clust.pi[[count0]] = pi.clust;
      log.like[[count0]] = loglik.iter;
      bic.clust[count0] = if(!junk.null.aware.bic){
        -2*tmp.loglik + (2*K.clust-1)*log(m);
      }else{
        -2*tmp.loglik + (2*K.clust - 1 - sum(junk.mixture) - sum(null.mixture))*log(m);
      }
      count0 = count0 +1;
    }
    bic.clust.mx[count.K] = max(bic.clust[(count0-(rand.num+2)):(count0-1)]);

    if((count.K > stop.bic.iter) & (i >= min.clust.search)){
      tmp.cond = TRUE;
      for(j in 1:stop.bic.iter){
        tmp.cond = tmp.cond & (bic.clust.mx[count.K-j+1] > bic.clust.mx[count.K-j]);
      }
      if(tmp.cond){
        break;
      }
    }
    count.K = count.K +1;
  }

  results = clust.prob(bic.clust, clust.mn, clust.pi, theta = theta, theta.sd = theta.se, junk.mixture = junk.mixture, df = df,
                       junk.mean = mu, junk.sd = sig, null.mixture = null.mixture, null.mean = mu.null, null.sd = sig.null,
                       obs.names = obs.names, clust.size.prior = clust.size.prior, prior = bic.prior);
  results.all = results[order(results$cluster),];
  results.best = best.clust(results)

  if(!is.null(cluster.membership)){
    variant.clusters = clust.inc.list(results, by.prob = cluster.membership$by.prob, bound = cluster.membership$bound);
  }
  if(!is.null(plot.results)){
    if(plot.results[[1]]=="best"){
      plot.1 = junk.clust.plot(results.best, sig = sig, mu = mu, mu.null = mu.null);
      plot.2 = two.stage.plot(results.best, bx, by, bxse, byse, obs.names);
    }else{
      tmp = plot.results[[2]];
      tmp.res = pr.clust(results, prob = tmp);
      plot.1 = junk.clust.plot(tmp.res, sig = sig, mu = mu, mu.null = mu.null, junk.mixture = junk.mixture, null.mixture = null.mixture);
      plot.2 = two.stage.plot(tmp.res, bx, by, bxse, byse, obs.names);
    }
    res = list(results = list(all = results.all, best = results.best),
               cluster.membership = variant.clusters, plots = list(two.stage = plot.2, split.plot = plot.1),
               log.likelihood = log.like, bic = bic.clust)
  }else{
    res = list(results = list(all = results.all, best = results.best),
               cluster.membership = variant.clusters,
               log.likelihood = log.like, bic = bic.clust)
  }
  if(trait.search){
    trt.search = pheno.search(variant.clusters = variant.clusters, p.value = trait.pvalue, r2 = proxy.r2, catalogue = catalogue, proxies = proxies, build = build);
    res = list(results = list(all = results.all, best = results.best),
               cluster.membership = variant.clusters, trait.search = trt.search,
               log.likelihood = log.like, bic = bic.clust);
  }
  return(res);
}
