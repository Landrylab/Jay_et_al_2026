args <- commandArgs(trailingOnly=TRUE)
out <- args[[1]]
conditional <- function(data, fields) {
  data <- data[data$carrier >= 0,,drop=FALSE]
  empty <- c(p=NA,or=NA,low=NA,high=NA,informative_strata=0,informative_samples=0)
  if (!nrow(data)) return(empty)
  strata <- interaction(data[,fields,drop=FALSE], drop=TRUE, lex.order=TRUE)
  tabs <- lapply(split(data,strata), function(d) table(factor(d$carrier,levels=c(1,0)),factor(d$phenotype,levels=c('R','S'))))
  tabs <- Filter(function(t) all(rowSums(t)>0) && all(colSums(t)>0), tabs)
  if (!length(tabs)) return(empty)
  a <- array(unlist(tabs),dim=c(2,2,length(tabs)))
  fit <- if(length(tabs)==1) fisher.test(a[,,1]) else mantelhaen.test(a,exact=TRUE,alternative='two.sided')
  c(p=fit$p.value,or=unname(fit$estimate),low=fit$conf.int[1],high=fit$conf.int[2],informative_strata=length(tabs),informative_samples=sum(a))
}
adjust <- function(p) {
  missing<-is.na(p);p[missing]<-1;q<-p.adjust(p,'BH');q[missing]<-NA;q
}
for (kind in c('gene','pool')) {
  calls <- read.csv(gzfile(file.path(out,paste0('sample_',kind,'_calls.csv.gz'))),stringsAsFactors=FALSE)
  result <- read.csv(file.path(out,paste0(kind,'_results.csv')),stringsAsFactors=FALSE)
  independent_p <- vapply(seq_len(nrow(result)), function(i) {
    tab <- matrix(as.numeric(unlist(result[i,c('R_carrier','R_noncarrier','S_carrier','S_noncarrier')])),2,byrow=TRUE)
    if(any(rowSums(tab)==0) || any(colSums(tab)==0)) return(NA_real_)
    fisher.test(tab)$p.value
  }, numeric(1))
  stopifnot(all(is.na(independent_p)==is.na(result$fisher_p)))
  stopifnot(max(abs(independent_p-result$fisher_p),na.rm=TRUE)<1e-10)
  stopifnot(max(abs(adjust(independent_p)-result$fisher_q_primary),na.rm=TRUE)<1e-10)
  key <- if(kind=='gene') 'gene_id' else 'pool'
  groups <- split(calls,paste(calls$drug,calls[[key]],sep='|'))
  for (i in seq_len(nrow(result))) {
    d <- groups[[paste(result$drug[i],result[[key]][i],sep='|')]]
    for(prefix in c('pca_cluster','pca_cluster_study')) {
      fields <- if(prefix=='pca_cluster') 'cluster_label' else c('cluster_label','Study_clean')
      vals <- conditional(d,fields)
      for(k in names(vals)) result[i,paste0(prefix,'_',k)]<-vals[k]
    }
  }
  for(prefix in c('pca_cluster','pca_cluster_study')) {
    ps<-paste0(prefix,'_p')
    result[[paste0(prefix,'_q_primary')]]<-adjust(result[[ps]])
    if(kind=='gene') {
      result[[paste0(prefix,'_q_original_edges')]]<-NA
      ix<-result$original_network_edge %in% c(TRUE,'True')
      result[ix,paste0(prefix,'_q_original_edges')]<-adjust(result[[ps]][ix])
    }
  }
  write.csv(result,file.path(out,paste0(kind,'_results.csv')),row.names=FALSE,na='')
}
