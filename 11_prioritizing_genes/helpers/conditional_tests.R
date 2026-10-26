args<-commandArgs(trailingOnly=TRUE);out<-args[1]
specs<-read.csv(file.path(out,'model_specifications.csv'),stringsAsFactors=FALSE)
samples<-read.csv(file.path(out,'model_sample_calls.csv'),stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(specs$model_id))
stopifnot(!anyDuplicated(samples[,c('model_id','sample')]))
stopifnot(all(samples$model_id %in% specs$model_id))
conditional<-function(d,fields) {
  empty<-list(p=NA_real_,or=NA_real_,low=NA_real_,high=NA_real_,informative_strata=0,informative_samples=0,informative_R=0,informative_S=0,status='not_estimable',diagnostic='No informative strata')
  if(!nrow(d))return(empty)
  key<-do.call(paste,c(d[,fields,drop=FALSE],sep='|'))
  groups<-split(d,key)
  tabs<-lapply(groups,function(z)table(factor(z$carrier,levels=c(1,0)),factor(z$phenotype,levels=c('R','S'))))
  tabs<-Filter(function(t)all(rowSums(t)>0)&&all(colSums(t)>0),tabs)
  if(!length(tabs))return(empty)
  a<-array(unlist(tabs),dim=c(2,2,length(tabs)))
  warnings<-character()
  fit<-tryCatch(withCallingHandlers(if(length(tabs)==1)fisher.test(a[,,1],alternative='two.sided') else mantelhaen.test(a,exact=TRUE,alternative='two.sided'),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(fit,'error')) {empty$diagnostic<-conditionMessage(fit);empty$status<-'error';return(empty)}
  est<-unname(fit$estimate)
  list(p=fit$p.value,or=est,low=fit$conf.int[1],high=fit$conf.int[2],informative_strata=length(tabs),informative_samples=sum(a),informative_R=sum(a[,1,]),informative_S=sum(a[,2,]),status=if(!is.finite(est)||est==0)'boundary_estimate' else 'estimable',diagnostic=paste(warnings,collapse=';'))
}
ans<-list()
for(i in seq_len(nrow(specs))) {
  s<-specs[i,,drop=FALSE];d<-samples[samples$model_id==s$model_id,,drop=FALSE]
  for(method in c('pca_only_same_samples','pca_joint')) {
    fields<-if(method=='pca_joint')c('cluster','partner_pattern') else 'cluster'
    ans[[length(ans)+1]]<-cbind(s,method=method,as.data.frame(conditional(d,fields),stringsAsFactors=FALSE))
  }
}
r<-do.call(rbind,ans);r$q_36_tests<-NA_real_
for(t in unique(r$threshold))for(m in unique(r$method)) {
  ix<-which(r$threshold==t&r$method==m);stopifnot(length(ix)==length(unique(specs$model_id[specs$threshold==t])))
  missing<-is.na(r$p[ix]);p<-r$p[ix];p[missing]<-1;q<-p.adjust(p,'BH');q[missing]<-NA_real_;r$q_36_tests[ix]<-q
}
# Preserve the exact threshold labels after the numerical grouping and tests.
r$threshold<-vapply(strsplit(r$model_id,'|',fixed=TRUE),tail,character(1),n=1)
write.csv(r,file.path(out,'conditional_results.csv'),row.names=FALSE,na='')
