args<-commandArgs(trailingOnly=TRUE)
out<-args[1]
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

m<-read.csv(file.path(out,'variant_metadata.csv'))
calls<-read.csv(file.path(out,'sample_variant_calls.csv.gz'))
groups<-split(calls,calls$key)
for(i in seq_len(nrow(m))) {
 d<-groups[[m$key[i]]]
 for(ph in c('R','S')) {
  z<-d[d$phenotype==ph,]
  m[i,paste0(ph,'_carrier')]<-sum(z$carrier==1)
  m[i,paste0(ph,'_called')]<-sum(z$carrier>=0)
 }
 vals<-conditional(d,'cluster_label')
 for(k in names(vals)) m[i,paste0('pca_',k)]<-vals[k]
}
m$pca_q_within_gene_drug<-NA_real_
for(fam in unique(m$family)) {
 ix<-which(m$family==fam);m$pca_q_within_gene_drug[ix]<-adjust(m$pca_p[ix])
}
write.csv(m,file.path(out,'variant_results.csv'),row.names=FALSE,na='')
