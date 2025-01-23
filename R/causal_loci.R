#' Identification of putative causal loci
#'
#' Identification of putative causal loci using KnockoffTrio's feature statistics
#' @param window The result window from function KnockoffTrio. If there are multiple result windows (e.g., when you analyze multiple regions in the genome), please use rbind to combine all the windows before running causal_loci.
#' @param M A positive integer for the number of knockoffs. The default is 10.
#' @param fdr A real number in a range of (0,1) indicating the target FDR level. The default is 0.1. Use 0.2 for a more lenient FDR control.
#' @return A list that contains the following elements for claiming significance using knockoff statistics. The result window also contains FBAT p-values and ACAT-combined p-values, which can be used for claiming significance in addition to knockoff statistics. If p-values are used, Bonferroni correction is usually necessary to adjust for multiple testing for controlling the family-wise error rate - see examples below. 
#' \describe{
#'   \item{window}{A data frame for an updated result window that includes an extra column for KnockoffTrio's Q-values. A locus with a Q-value <= the target FDR level, i.e., window$q<=fdr, is considered as putative causal at the target FDR.}
#'   \item{thr.w}{A positive real number indicating the significance threshold for KnockoffTrio's feature statistics. A locus with a feature statistic >= thr.w, i.e., window$w>=thr.w is considered as putative causal at the target FDR. The loci selected by window$w>=thr.w are equivalent to those by window$q<=fdr. No loci are selected at the target FDR level if thr.w=Inf.}
#' }
#' @importFrom stats as.dist cutree hclust median pcauchy pnorm princomp
#' @export
#' @examples
#' data(KnockoffTrio.example)
#' knockoff<-create_knockoff(trio.hap=KnockoffTrio.example$trio.hap,
#'           duo.hap=KnockoffTrio.example$duo.hap, pos=KnockoffTrio.example$pos, M=10)
#' window<-KnockoffTrio(trio=KnockoffTrio.example$trio, trio.ko=knockoff$trio.ko,
#'         duo=knockoff$duo, duo.ko=knockoff$duo.ko, pos=KnockoffTrio.example$pos)
#' 
#' #Identification of significant loci using KnockoffTrio's feature statistics (W or Q) 
#' #at a target FDR
#' target_fdr<-0.1
#' result<-causal_loci(window,M=10,fdr=target_fdr)
#' sig_loci_by_w_index<-which(result$window$w>=result$thr.w)
#' sig_loci_by_q_index<-which(result$window$q<=target_fdr)
#' 
#' #Identification of significant loci using FBAT p-values with Bonferroni correction
#' #for controlling the family-wise error rate at 0.05
#' sig_loci_by_p_fbat_index<-which(window$p.burden<0.05/nrow(window))
#' 
#' #Identification of significant loci using ACAT p-values with Bonferroni correction
#' #for controlling the family-wise error rate at 0.05
#' sig_loci_by_p_acat_index<-which(window$p<0.05/nrow(window))
causal_loci<-function(window,M=10,fdr=0.1){
  q<-MK.q.byStat(window$kappa,window$tau,M=M)
  window<-cbind(window,q)
  thr.w<-MK.threshold.byStat(window$kappa,window$tau,M=M,fdr=fdr)
  out<-list()
  out$window<-window
  out$thr.w<-thr.w
  return(out)
}
MK.q.byStat<-function (kappa,tau,M,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau))
  for(i in 1:length(b)){
    q[b[i]]<-min(ratio[i:min(length(b),Rej.Bound)])*c_0[i]+1-c_0[i]
    if(i>Rej.Bound){break}
  }
  return(q)
}
MK.threshold.byStat<-function (kappa,tau,M,fdr = 0.1,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  ok<-which(ratio<=fdr)
  if(length(ok)>0){
    return(tau[b][ok[length(ok)]])
  }else{return(Inf)}
}

