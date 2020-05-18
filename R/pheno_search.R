## search phenoscanner
pheno_search <- function(variant_clusters, p_value = 1e-5, r2 = 0.8,
                         catalogue = "GWAS", proxies = "None", build = 37) {
  tmp_search <- names(variant_clusters)
  trait_search <- data.frame(matrix(vector(), 0, 5, dimnames = list(c(),
        c("cluster", "trait", "associated.snps", "total.snps", "p.cut.off"))),
        stringsAsFactors = F)
  trt_count <- 1
  for (i_clst in seq_len(length(tmp_search))) {
    tmp_nam <- strsplit(tmp_search[i_clst], split = "cluster_")[[1]]
    tmp <- !tmp_nam[2] == "Null" & !tmp_nam[2] == "Junk"
    if (tmp) {
      tmp_snps <- variant_clusters[[i_clst]][1]
      phenodata <- MendelianRandomization::phenoscanner(snpquery =
                      tmp_snps[[1]], pvalue = p_value, r2 = r2,
                      catalogue = catalogue, proxies = proxies, build = build)
      traits <- unique(phenodata$results$trait)
      for (i_trts in seq_len(length(traits))) {
        pos <- unique(phenodata$results$snp[phenodata$results$trait ==
                                              traits[i_trts]])
        trait_search[trt_count, ] <- c(tmp_search[i_clst], traits[i_trts],
                                       length(pos), length(tmp_snps[[1]]),
                                       p_value)
        trt_count <- trt_count + 1
      }
    }
  }
  return(trait_search)
}
