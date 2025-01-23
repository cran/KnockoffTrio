#' Meta-analysis for KnockoffTrio
#'
#' Meta-analysis for KnockoffTrio
#'
#' @param window A list of windows for the analysis results from different cohorts/studies.
#' @param n A positive integer vector for the number of families in each cohort/study. For weighted meta-analysis, a study's weight is based on the number of families. The default is NA for unweighted meta-analysis.
#' @param M A positive integer for the number of knockoffs. The default is 10.
#' @return A data frame for the meta-analysis results.
#' @importFrom stats as.dist cutree hclust median pcauchy pnorm princomp
#' @export
#' @examples
#' data(KnockoffTrio.example)
#' knockoff<-create_knockoff(trio.hap=KnockoffTrio.example$trio.hap,
#'           duo.hap=KnockoffTrio.example$duo.hap, pos=KnockoffTrio.example$pos, M=10)
#' window<-KnockoffTrio(trio=KnockoffTrio.example$trio, trio.ko=knockoff$trio.ko,
#'         duo=knockoff$duo, duo.ko=knockoff$duo.ko, pos=KnockoffTrio.example$pos)
#' window.list<-list(window,window)
#' window.meta<-meta_analysis(window.list,M=10)
#' result<-causal_loci(window.meta,M=10,fdr=0.1)
meta_analysis<-function(window,n=NA,M=10){
  position<-c()
  nwin<-length(window)
  if (is.na(n)) n<-rep(1,nwin)
  w<-sqrt(n)
  win<-vector(mode = "list", length = nwin)
  for (i in 1:nwin) {
    win[[i]]<-paste0(window[[i]][,"start"],window[[i]][,"end"])
    position<-c(position,win[[i]])
  }
  position<-sort(unique(position))
  nvar<-length(position)
  ind<-vector(mode = "list", length = nwin)
  for (i in 1:nwin) ind[[i]]<-match(position,win[[i]])
  #out<-array(NA,dim=c(nvar,ncol(window[[1]])))
  out<-window[[1]][-(1:nrow(window[[1]])),]
  #colnames(out)<-colnames(window[[1]])
  coln<-paste0("z_",1:M)
  colnp<-paste0("p_",1:M)
  zcoln<-c("z",coln)
  pcoln<-c("p",colnp)
  for (pos in 1:nvar){
    nmeta<-c()
    for (i in 1:nwin){
      if (!is.na(ind[[i]][pos])){
        #meta[[pos]]<-c(meta[[pos]],i)
        if (length(nmeta)==0){
          out[pos,]<-window[[i]][ind[[i]][pos],]
          out[pos,zcoln]<-out[pos,zcoln]*w[i]
        }else{
          out[pos,zcoln]<-out[pos,zcoln]+window[[i]][ind[[i]][pos],zcoln]*w[i]
        }
        nmeta<-c(nmeta,i)
      }
    }
    denom<-sqrt(sum(n[nmeta]))
    out[pos,zcoln]<-out[pos,zcoln]/denom
    out[pos,pcoln]<-2*pnorm(-abs(as.matrix(out[pos,zcoln])))
  }
  tmp<-calculate_w_kappatau(out[,"p"],t(out[,colnp]))
  out[,"w"]<-tmp$w
  out[,"kappa"]<-tmp$kappatau[,1]
  out[,"tau"]<-tmp$kappatau[,2]
  out<-as.data.frame(out)
  return(out)
}
calculate_w_kappatau<-function(q1,q2){
  out<-list()
  t1<--log10(q1)
  t2<--log10(q2)
  t2_med<-apply(t2,2,median)
  t2_max<-apply(t2,2,max)
  out$w<-(t1-t2_med)*(t1>=t2_max)
  out$w.raw<-t1-t2_med
  out$kappatau<-MK.statistic(t1,t(t2),method="median")
  #out$q<-MK.q.byStat(out$kappatau[,1],out$kappatau[,2],M=M)
  return(out)
}
MK.statistic<-function (T_0,T_k,method='median'){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0
  kappa<-apply(T.temp,1,which.max)-1
  if(method=='median'){
    Get.OtherMedian<-function(x){median(x[-which.max(x)])}
    tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  }
  return(cbind(kappa,tau))
}
