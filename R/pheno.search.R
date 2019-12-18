## search phenoscanner
pheno.search = function(variant.clusters, p.value = 1e-5, r2 = 0.8, catalogue = "GWAS", proxies = "None", build = 37){
  tmp.search = names(variant.clusters);
  trait.search = data.frame(matrix(vector(), 0,5, dimnames=list(c(), c("cluster", "trait", "associated.snps", "total.snps", "p.cut.off"))), stringsAsFactors=F)
  trt.count = 1
  for(i.clst in 1:length(tmp.search)){
    tmp.nam = strsplit(tmp.search[i.clst], split = "cluster_")[[1]];
    tmp = !tmp.nam[2] == "Null"  & !tmp.nam[2] == "Junk";
    if(tmp){
      tmp.snps = variant.clusters[[i.clst]][1];
      phenodata = MendelianRandomization::phenoscanner(snpquery = tmp.snps[[1]], pvalue=p.value, r2 = r2, catalogue = catalogue, proxies = proxies, build = build);
      traits = unique(phenodata$results$trait);
      for(i.trts in 1:length(traits)){
        pos = unique(phenodata$results$snp[phenodata$results$trait==traits[i.trts]]);
        trait.search[trt.count, ] = c(tmp.search[i.clst], traits[i.trts], length(pos), length(tmp.snps[[1]]), p.value);
        trt.count = trt.count + 1;
      }
    }
  }
  return(trait.search)
}
