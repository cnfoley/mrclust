## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=5, 
  fig.height= 4
)

## ----setup--------------------------------------------------------------------
library(mrclust)

## ---- echo=F------------------------------------------------------------------
options(warn=-1)

## ---- eval=F------------------------------------------------------------------
#  install.packages("devtools", repos='http://cran.us.r-project.org')
#  library(devtools)
#  install_github("cnfoley/mrclust", build_vignettes = TRUE)

## -----------------------------------------------------------------------------
library(mrclust)
sbp.cad <- mrclust::SBP_CAD
head(sbp.cad)

## -----------------------------------------------------------------------------
# dbp.cad <- mrclust::DBP_CAD
# pp.cad <- mrclust::PP_CAD

## -----------------------------------------------------------------------------
?SBP_CAD

## -----------------------------------------------------------------------------
##############################################################
# estimates of the G-X associations and their standard errors
##############################################################

bx = sbp.cad$bx;
bxse = sbp.cad$bxse;

##############################################################
# estimates of the G-Y associations and their standard errors
##############################################################

by = sbp.cad$by;
byse = sbp.cad$byse;

###########################################
# ratio-estimates and their standard errors
###########################################

ratio.est = by/bx;
ratio.est.se = byse/abs(bx);

#####################
# Names for the snps
#####################
snp.names = sbp.cad$chr.pos;

## -----------------------------------------------------------------------------
############################
########### SBP-CAD analysis
############################

res_em = mr_clust_em(theta=ratio.est, theta.se=ratio.est.se, bx=bx, by=by, bxse=bxse, byse=byse,
                       obs.names = snp.names);


## -----------------------------------------------------------------------------
names(res_em);

## -----------------------------------------------------------------------------
#############################
# results contain two tables:
#############################

names(res_em$results);

####################################################################################################
# first few rows of table containing summaries for the genetic variants allocated to "all" clusters 
####################################################################################################

head(res_em$results$all);

#######################################################################################################
# first few rows of table containing summaries for the genetic variants allocated to the "best" cluster
#######################################################################################################

head(res_em$results$best);


## -----------------------------------------------------------------------------
######################
# Identified clusters
######################

clusters = unique(res_em$results$best$cluster.class); # same answer if we type unique(res_em$results$all$cluster.class)
clusters[order(clusters)];

######################
# number of clusters K
######################

K = length(clusters);
K;

## -----------------------------------------------------------------------------
######################
#  clusters
######################

names(res_em$cluster.membership)

###########################################################################
# Variants assinged to the first few deciles of the allocation probability 
###########################################################################

head(res_em$cluster.membership$cluster_1);

## -----------------------------------------------------------------------------
#########################################################
# Variants assinged to cluster one with probability >0.9 
#########################################################

res_em$cluster.membership$cluster_1$pr_range_1_to_0.9;

## -----------------------------------------------------------------------------
###################################
# Scatter plot of the fitted model  
###################################

plot.sbp.best = res_em$plots$two.stage + ggplot2::xlim(0, max(abs(bx)+2*bxse)) + ggplot2::xlab("Genetic association with SBP") + ggplot2::ylab("Genetic association with CAD") + ggplot2::ggtitle("");

##############
# plot output
##############

plot.sbp.best

## -----------------------------------------------------------------------------
#################################################################################################
# Variants allocated to clusters with: (i) at least 4 variants and; (ii) allocation probability >=0.8 
#################################################################################################

res80 = mrclust::pr.clust(dta = res_em$results$best, prob = 0.8, min.obs =  4);

## -----------------------------------------------------------------------------
##############################################################
# Update summary effects using the prioritised set of variants  
##############################################################

      keep80 = which(snp.names %in% res80$observation);
      bx80   = bx[keep80];
      bxse80 = bxse[keep80];
      by80   = by[keep80];
      byse80 = byse[keep80];
      snp.names80 = snp.names[keep80];

plot.sbp.pr80 = two.stage.plot(res = res80, bx = bx80, by = by80, bxse = bxse80, byse = byse80, obs.names = snp.names80) + ggplot2::xlim(0, max(abs(bx80)+2*bxse80)) + ggplot2::xlab("Genetic association with SBP") + ggplot2::ylab("Genetic association with CAD") + ggplot2::ggtitle("");
# plot result
plot.sbp.pr80;

## -----------------------------------------------------------------------------
###################################################
# Re-run analysis with trait search option == TRUE
###################################################

res_em = mr_clust_em(theta=ratio.est, theta.se=ratio.est.se, bx=bx, by=by, bxse=bxse, byse=byse,
                       obs.names = snp.names, max.iter = 5e3, tol = 1e-5,
                       junk.mean = 0,
                       results.list = list("all", "best"),
                       cluster.membership = list(by.prob = 0.25, bound = 0.75),
                       plot.results = list("best", min.pr = 0.5),
                       trait.search = TRUE, trait.pvalue = 1e-5, proxy.r2 = 0.8,
                       catalogue = "GWAS", proxies = "None", build = 37);


## -----------------------------------------------------------------------------
##############################
# Access trait search results
##############################

res.search = res_em$trait.search;

####################
# variables returned
####################

names(res.search);

## -----------------------------------------------------------------------------
############################################################
# Setting
# "cluster.membership = list(by.prob = 0.25, bound = 0.75)"
# in the argument of "mr_clust_em"
# returns:
############################################################

head(res_em$cluster.membership);

## -----------------------------------------------------------------------------
###############################################
# Count the number of variants in each cluster 
###############################################

tab.clust = table(res_em$results$best$cluster.class)
tab.clust;

##########################################################
# Identify the cluster with the largest number of variants 
##########################################################

mx.clust = which.max(tab.clust[!names(tab.clust) %in% c("Null", "Junk")]);
clst.trt.srch = paste0("cluster_",mx.clust) 


## -----------------------------------------------------------------------------
##############################
# Half the total number of snp
##############################

diff.snps.tol = floor(tab.clust[mx.clust]/2);

##############################################
# Difference between total and associated snps
##############################################

diff.snps = as.numeric(res_em$trait.search$total.snps)-as.numeric(res_em$trait.search$associated.snps);


## -----------------------------------------------------------------------------
head(res_em$trait.search[which(diff.snps<=diff.snps.tol & res_em$trait.search$cluster==clst.trt.srch),]);

