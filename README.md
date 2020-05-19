# mrclust
R package for investigating clustered heterogeneity in Mendelian randomization
(MR) analyses.

Genetic variants which recover similar estimates of the causal effect of the
risk-factor on the outcome - i.e. their ratio-estimates are similar in
direction, magnitude and precision - can form distinct clusters in MR analyses.
We call this `clustered heterogeneity'. Each cluster might represent a distinct
pathway with which the risk-factor is related to the outcome and hence clustered
heterogeneity is interesting to investigate as the identity of the genetic
variants in the clusters may reveal information about the risk factor and how it
relates to the outcome. 

## Functions
* mr_clust_em - performs expectation-maximisation (EM) based model fitting of
the MR-Clust mixture model.

## Installation
1. install.packages("devtools")
2. library(devtools)
3. install_github("cnfoley/mrclust", build_vignettes = TRUE)
4. library(mrclust)
5. browseVignettes("mrclust")


## Example
\# Regression coefficients and standard errors from systolic blood pressure and
coronary artery disease studies.
* sbp_cad = mrclust::SBP_CAD
* bx = sbp_cad$bx
* bxse = sbp_cad$bxse
* by = sbp_cad$by
* byse = sbp_cad$byse
* ratio_est = by/bx
* ratio_est_se = byse/abs(bx)


\# SNP IDs  
* snp_names = sbp_cad$chr.pos;

### Clustering analysis
* res_em = mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx,
by = by, bxse = bxse, byse = byse, obs_names = snp_names)
### Table of variant-cluster allocations
* head(res_em$results$best)

### A cluster annotated scatter-plot of the G-X and G-Y associations
* plot.sbp.best = res_em$plots$two.stage +
ggplot2::xlim(0, max(abs(bx)+2*bxse)) +
ggplot2::xlab("Genetic association with SBP") +
ggplot2::ylab("Genetic association with CAD") +
ggplot2::ggtitle("")
* plot.sbp.best
