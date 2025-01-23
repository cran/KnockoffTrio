#' Calculate KnockoffTrio's feature statistics
#'
#' Calculate KnockoffTrio's feature statistics and FBAT statistics using original and knockoff genotype data.
#'
#' @param trio A 3n*p matrix for the trio genotype data, in which n is the number of trios and p is the number of variants. Each trio must consist of father, mother, and offspring (in this order). The genotypes must be coded as 0, 1, or 2. Missing genotypes are not allowed.
#' @param trio.ko A 3n*p*M array for the knockoff trio genotype data created by function create_knockoff. M is the number of knockoffs.
#' @param duo A 3m*p matrix for the duo genotype data created by function create_knockoff, in which m is the number of duos and p is the number of variants. Please do not use the original 2m*p duo genotype matrix.
#' @param duo.ko A 3m*p*M array for the knockoff duo genotype data created by function create_knockoff. M is the number of knockoffs.
#' @param pos A numeric vector of length p for the position of p variants.
#' @param start An integer for the first position of sliding windows. If NULL, start=min(pos). Only used if you would like to use the same starting position for different cohorts/analyses.
#' @param end An integer for the last position of sliding windows. If NULL, end=max(pos). Only used if you would like to use the same ending position for different cohorts/analyses.
#' @param size A numeric vector for the size(s) of sliding windows when scanning the genome.
#' @param p_value_only A logical value indicating whether to perform the knockoff analysis. When p_value_only is TRUE, only the ACAT-combined p-values are to be calculated for each window. When p_value_only is FALSE, trio.ko or duo.ko is required and KnockoffTrio's feature statistics are to be calculated for each window in addition to the p-values.
#' @param adjust_for_cov A logical value indicating whether to adjust for covariates. When adjust_for_cov is TRUE, y is required.
#' @param y A numeric vector of length n for the residual Y-Y_hat. Y_hat is the predicted value from the regression model in which the quantitative trait Y is regressed on the covariates. If Y is dichotomous, you may treat Y as quantitative when applying the regression model.
#' @param chr A character for the name of the chromosome, e.g., "1", "2", ..., "22", and "X".
#' @param xchr A logical value indicating whether the analysis is for the X chromosome. When xchr is TRUE, the analysis is for the X chromosome and sex is required. When xchr is FALSE, the analysis is for the autosomes. The default if FALSE.
#' @param sex A numeric vector of length n for the sex of offspring. 0s indicate females and 1s indicate males. Sex is required when xchr is TRUE.
#' @return A data frame for analysis results from KnockoffTrio and FBAT. The data frame contains the following columns if p_value_only is FALSE:
#' \describe{
#'   \item{chr}{The chromosome number.}
#'   \item{start, end}{The start and end position of a window.}
#'   \item{actual_start, actual_end}{The position of the first and last variant in a window.}
#'   \item{n}{The number of variants in a window.}
#'   \item{dir}{The direction of effect of the most significant variant in a window.}
#'   \item{w}{The W knockoff feature statistic for a window. Please use function causal_loci to obtain the significance threshold for w at target FDRs.}
#'   \item{p}{The ACAT-combined p-value for a window. If a window contains multiple variants (i.e., n>1), ACAT combines FBAT p-values for each variant and a burden FBAT p-value for all variants in the window. If a window contains only one variant (i.e., n=1), the ACAT-combined p-value is equivalent to the FBAT p-value for this variant.}
#'   \item{z}{The FBAT z-score for a window. If a window contains multiple variants (i.e., n>1), z is the burden FBAT z-score for all variants in the window. If a window contains only one variant (i.e., n=1), z is the FBAT z-score for this variant.}
#'   \item{p.burden}{The FBAT p-value for a window. If a window contains multiple variants (i.e., n>1), p.burden is the burden FBAT p-value for all variants in the window. If a window contains only one variant (i.e., n=1), p.burden is the FBAT p-value for this variant.}
#'   \item{kappa, tau}{The two columns are used by function causal_loci for knockoff inference.}
#'   \item{p_1, ..., p_M}{The ACAT-combined p-values for M knockoffs.}
#'   \item{z_1, ..., z_M}{The FBAT z-scores for M knockoffs.}
#' }
#' @importFrom stats as.dist cutree hclust median pcauchy pnorm princomp
#' @export
#' @examples
#' data(KnockoffTrio.example)
#' knockoff<-create_knockoff(trio.hap=KnockoffTrio.example$trio.hap,
#'           duo.hap=KnockoffTrio.example$duo.hap, pos=KnockoffTrio.example$pos, M=10)
#'
#' #Analysis for both trios and duos
#' window<-KnockoffTrio(trio=KnockoffTrio.example$trio, trio.ko=knockoff$trio.ko,
#'         duo=knockoff$duo, duo.ko=knockoff$duo.ko, pos=KnockoffTrio.example$pos)
#'
#' #Analysis for trios only
#' window<-KnockoffTrio(trio=KnockoffTrio.example$trio, trio.ko=knockoff$trio.ko,
#'         duo=NULL, duo.ko=NULL, pos=KnockoffTrio.example$pos)
#'
#' #Analysis for duos only
#' window<-KnockoffTrio(trio=NULL, trio.ko=NULL,
#'         duo=knockoff$duo, duo.ko=knockoff$duo.ko, pos=KnockoffTrio.example$pos)
KnockoffTrio<-function(trio, trio.ko=NULL, duo=NULL, duo.ko=NULL, pos, start=NULL, end=NULL, size=c(1,1000,5000,10000,20000,50000), p_value_only=FALSE, adjust_for_cov=FALSE, y=NULL, chr="1", xchr=FALSE, sex=NULL){
  if (!is.null(trio)) {
    if (nrow(trio) %% 3!=0) stop("The number of rows of the trio genotype matrix must be a multiple of three.")
  }
  if (!is.null(duo)) {
    if (nrow(duo) %% 3!=0) stop("The number of rows of the duo genotype matrix must be a multiple of three due to FBAT. Please use the duo genotype matrix generated from the create_knockoff function.")
  }
  if (adjust_for_cov & is.null(y)) stop("When adjust_for_cov is TRUE, y cannot be NULL.")
  if (!is.null(trio) & !is.null(duo)){
    dat<-rbind(trio, duo)
    if (!p_value_only){
      if (is.null(trio.ko) | is.null(duo.ko)) stop("The trio and duo knockoff data cannot be NULL if knockoff analysis is to be performed on both trio and duo data. Please use the create_knockoff function to generate knockoff data.")
      if (dim(trio.ko)[3]!=dim(duo.ko)[3]) stop("The trio and duo knockoff data must have the same number of knockoffs.")
      M<-dim(trio.ko)[3]
      dat.ko<-array(dim=c(dim(dat), M))
      for (i in 1:M) dat.ko[,,i]<-rbind(trio.ko[,,i,drop=F], duo.ko[,,i,drop=F])
    } else dat.ko<-NULL
    rm(trio)
    rm(trio.ko)
    rm(duo)
    rm(duo.ko)
  } else if (!is.null(trio)){
    dat<-trio
    dat.ko<-trio.ko
    rm(trio)
    rm(trio.ko)
  } else if (!is.null(duo)){
    dat<-duo
    dat.ko<-duo.ko
    rm(duo)
    rm(duo.ko)
  }
  #if (nrow(dat) %% 3!=0) stop("The number of rows of the original trio matrix must be a multiple of three.")
  if (xchr & is.null(sex)) stop("Gender information is required if KnockoffTrio is applied to the X chromosome.")
  if (is.null(start)) start<-min(pos)
  if (is.null(end))  end<-max(pos)

  n<-nrow(dat)/3 #number of trios
  nsnp<-ncol(dat)
  out_fbat<-fbat(dat,adjust_for_cov,y,xchr=xchr,sex=sex)
  z<-out_fbat$additive
  p1<-2*pnorm(-abs(z)) #p-values for SNP-based FBAT in original trios
  p1_sign<-sign(z)

  index_dad<-seq(1,(3*n),3) #index for parent1
  index_mom<-seq(2,(3*n),3) #index for parent2
  index_off<-seq(3,(3*n),3) #index for offspring
  mac<-apply(dat[-index_off,],2,sum) #minor allele count for parents
  if (xchr==FALSE) maf<-mac/(4*n) else maf<-(apply(dat[index_dad,],2,sum)/n+apply(dat[index_mom,],2,sum)/(2*n))/2

  #Prepare Z matrix for set-based fbat
  weight<-1/(sqrt(n*maf*(1-maf)))
  weight[which(weight==Inf)]<-0
  #original
  Z<-sweep(out_fbat$Z, 2, weight, '*')
  Zsq<-t(Z)%*%Z
  Zsum<-apply(Z,2,sum)

  #Prepare windows
  window<-data.frame()
  index<-which(maf>=0) #index<-which(maf>=0.01)
  for (i in 1:length(size)){
    if (size[i]==1) {
      if (length(index)>0) window<-data.frame(start=pos[index],end=pos[index],ind3=index,ind4=index,n=rep(1,length(index)))
    } else window<-rbind(window,getslidingwindow(start,end,size[i],pos))
  }
  window<-window[which(window$n>0),]
  window<-window[!duplicated(window[,c(3,4)]),]
  tmp<-which(window$n==1 & !(window$ind3 %in% index))
  if (length(tmp)>0) window<-window[-tmp,]
  nwindow<-nrow(window)
  q1<-p.fbat<-z1<-rep(NA,nwindow) #window p-values for original trios q1.cv<-q1.rv<-q1.pa<-q1.ma<-
  #q1.burden.cv<-q1.burden.rv<-q1.burden.urv<-rep(NA,nwindow) #p.denovo
  #n.denovo<-rep(0,nwindow)
  direction<-rep(NA,nwindow) #direction.pa<-direction.ma<-

  if (!p_value_only){
    #The following three warnings only work when xchr==FALSE
    #if (dim(dat.ko)[1] %% 3!=0) stop("The number of rows of the knockoff trio matrix must be a multiple of three.")
    #if (nrow(dat)!=dim(dat.ko)[1]) stop("The number of rows of the two matrices must be equal.")
    #if (ncol(dat)!=dim(dat.ko)[2]) stop("The number of columns of the two matrices must be equal.")
    M<-dim(dat.ko)[3]
    out_fbat<-fbat(dat,adjust_for_cov,y,dosage=TRUE,dat.ko,xchr,sex) #SNP-based FBAT statistics for knockoff trios
    z.ko<-out_fbat$additive
    Z2<-out_fbat$Z
    p2<-2*pnorm(-abs(z.ko)) #p-values for SNP-based FBAT in knockoff trios, dim=c(nsnp,M)

    Z2sq<-array(dim=c(nsnp,nsnp,M))
    Z2sum<-array(dim=c(nsnp,M))
    if (dim(Z2)[2]>1){
      for (m in 1:M){
        Z2[,,m]<-sweep(Z2[,,m], 2, weight, '*')
        Z2sq[,,m]<-t(Z2[,,m])%*%Z2[,,m]
        Z2sum[,m]<-apply(Z2[,,m],2,sum)
      }
    } else if (dim(Z2)[2]==1){
      for (m in 1:M){
        Z2[,,m]<-Z2[,,m]*weight
        Z2sq[,,m]<-t(Z2[,,m])%*%Z2[,,m]
        Z2sum[,m]<-sum(Z2[,,m])
      }
    }
    q2<-z2<-array(dim=c(M,nwindow))
  }

  for (i in 1:nwindow) {
    ind<-window$ind3[i]:window$ind4[i]

    if (length(ind)==1){
      if (maf[ind]>=0.01){
        q1[i]<-p.fbat[i]<-p1[ind] #q1.cv[i]
        z1[i]<-z[ind]
        direction[i]<-p1_sign[ind]
      }
    }else{
      z1[i]<-fbat_set(Zsum,Zsq,ind)
      p.fbat[i]<-2*pnorm(-abs(z1[i]))

      #single-variant for all
      ind.single<-ind #[mac[ind]>=5]
      p.single<-p1[ind.single]
      q1[i]<-ACAT(c(p.single,p.fbat[i]))

      #direction of q1[i]
      stat<-c(z1[i],z[ind.single])
      tmp<-which.max(abs(stat))
      if (length(tmp)>0) direction[i]<-sign(stat[tmp])

    }

    if (!p_value_only){
      if (length(ind)==1){
        if (maf[ind]>=0.01){
          q2[,i]<-p2[ind,] #q2.cv[,i]<-
          z2[,i]<-z.ko[ind,]

        }
      }else{
        for (m in 1:M){
          p.single.ko<-p2[ind.single,m]

          z2[m,i]<-fbat_set(Z2sum[,m],Z2sq[,,m],ind)
          p.burden.ko<-2*pnorm(-abs(z2[m,i]))
          q2[m,i]<-ACAT(c(p.single.ko,p.burden.ko))
        }
      }
    }
  }

  window$ind3<-pos[window$ind3]
  window$ind4<-pos[window$ind4]
  colnames(window)[3]<-"actual_start"
  colnames(window)[4]<-"actual_end"
  if (!p_value_only){
    out<-calculate_w_kappatau(q1,q2)
    w<-out$w
    kappatau<-out$kappatau

    rownames(q2)<-paste0("p_",1:M)
    q2<-t(q2)
    rownames(z2)<-paste0("z_",1:M)
    z2<-t(z2)
    #p: acat p, z: z score for single-variant or set FBAT (depends on if a window contains 1 or more variants), p.burden: p value for single-variant or set FBAT
    window<-cbind(chr,window,dir=direction,w,p=q1,z=z1,p.burden=p.fbat,kappatau,q2,z2)
  } else window<-cbind(chr,window,dir=direction,p=q1,z=z1,p.burden=p.fbat)

  return(window)
}
fbat<-function(dat,adjust_for_cov=FALSE,y=NULL,dosage=FALSE,dat1=NULL,xchr=FALSE,sex=NULL){
  #if (xchr & is.null(sex)) stop("Gender information is required if FBAT is applied to the X chromosome")
  if (is.null(ncol(dat))) dat<-as.matrix(dat)
  n<-nrow(dat)/3
  index_dad<-seq(1,(3*n),3) #index for parent1
  index_mom<-seq(2,(3*n),3) #index for parent2
  index_off<-seq(3,(3*n),3) #index for offspring
  if (!xchr){
    if (dosage==FALSE) Z<-dat[index_off,]-(dat[index_dad,]+dat[index_mom,])/2 else Z<-dat1[index_off,,,drop=FALSE]-(dat1[index_dad,,,drop=FALSE]+dat1[index_mom,,,drop=FALSE])/2
  } else Z<-xcontribution(dat,dat1,dosage,sex)
  # tat<-pcontribution(dat,dat1,dosage,xchr,sex)
  # Z.pa<-tat$Z.pa
  # Z.ma<-tat$Z.ma
  if (adjust_for_cov){
    Z<-sweep(Z,1,y,'*')
    # Z.pa<-sweep(Z.pa,1,y,'*')
    # Z.ma<-sweep(Z.ma,1,y,'*')
  }
  additive<-colSums(Z)/sqrt(colSums(Z^2)) #if dosage=T, dim=c(nsnp,M)
  # paternal<-colSums(Z.pa)/sqrt(colSums(Z.pa^2)) #if dosage=T, dim=c(nsnp,M)
  # maternal<-colSums(Z.ma)/sqrt(colSums(Z.ma^2)) #if dosage=T, dim=c(nsnp,M)
  additive[is.na(additive)]<-0
  # paternal[is.na(paternal)]<-0
  # maternal[is.na(maternal)]<-0
  out<-list()
  out$additive<-additive
  # out$paternal<-paternal
  # out$maternal<-maternal
  out$Z<-Z
  # out$Z.pa<-Z.pa
  # out$Z.ma<-Z.ma
  return (out)
}
fbat_set<-function(W,V,ind){
  if (length(ind)>0){
    s1<-sum(W[ind])
    s2<-sum(V[ind,ind])
    if (s2==0) s2<-1
    return (unname(s1/sqrt(s2)))
  } else return(NA)
}
pcontribution <- function(dat,dat1=NULL,dosage=FALSE,xchr=FALSE,sex=NULL){
  n.row <- nrow(dat)
  dad <- dat[seq.int(1, n.row, 3),, drop=FALSE]
  mom <- dat[seq.int(2, n.row, 3),, drop=FALSE]
  kid <- dat[seq.int(3, n.row, 3),, drop=FALSE]
  if (dosage==FALSE){
    Z.pa<-Z.ma<-array(0,dim=c(n.row/3,ncol(dat)))
    het <- (mom == 1)
    hethom <- het & (dad == 2)
    # ind212<-hethom & (kid == 2)
    # ind211<-hethom & (kid == 1)
    # n212 <- colSums(ind212, na.rm=TRUE)
    # n211 <- colSums(ind211, na.rm=TRUE)
    Z.ma[hethom & (kid == 2)]<-1
    Z.ma[hethom & (kid == 1)]<-(-1)
    hethom <- het & (dad == 0)
    # ind011<-hethom & (kid == 1)
    # ind010<-hethom & (kid == 0)
    # n011 <- colSums(ind011, na.rm=TRUE)
    # n010 <- colSums(ind010, na.rm=TRUE)
    Z.ma[hethom & (kid == 1)]<-1
    Z.ma[hethom & (kid == 0)]<-(-1)
    het <- (dad == 1)
    hethom <- het & (mom == 2)
    # ind122<-hethom & (kid == 2)
    # ind121<-hethom & (kid == 1)
    # n122 <- colSums(ind122, na.rm=TRUE)
    # n121 <- colSums(ind121, na.rm=TRUE)
    Z.pa[hethom & (kid == 2)]<-1
    Z.pa[hethom & (kid == 1)]<-(-1)
    hethom <- het & (mom == 0)
    # ind101<-hethom & (kid == 1)
    # ind100<-hethom & (kid == 0)
    # n101 <- colSums(ind101, na.rm=TRUE)
    # n100 <- colSums(ind100, na.rm=TRUE)
    Z.pa[hethom & (kid == 1)]<-1
    Z.pa[hethom & (kid == 0)]<-(-1)

    het <- (mom == 1) & (dad == 1)
    ind112<-het & (kid == 2)
    ind110<-het & (kid == 0)
    # n112 <- colSums(ind112, na.rm=TRUE)
    # n110 <- colSums(ind110, na.rm=TRUE)
    Z.pa[ind112]<-Z.ma[ind112]<-1
    Z.pa[ind110]<-Z.ma[ind110]<-(-1)
    #mat<-cbind(n100,n110,n121,n101,n112,n122,n010,n110,n211,n011,n112,n212)
    #colnames(mat)<-c("p-1","p-1","p-1","p+1","p+1","p+1","m-1","m-1","m-1","m+1","m+1","m+1") #contributions

    if (xchr){
      ind <- (dad == 1) & (mom == 1) & (kid==1)
      ind0<-sweep(ind,1,sex==0,'&')
      Z.pa[ind0]<-1
      Z.ma[ind0]<-(-1)
      ind0<-sweep(ind,1,sex==1,'&')
      Z.pa[ind0]<-(-1)
      Z.ma[ind0]<-1
    }

    out<-list()
    out$Z.pa<-Z.pa
    out$Z.ma<-Z.ma
    return(out)
  }else{
    dims<-dim(dat1)
    dims[1]<-dims[1]/3
    Z.pa<-Z.ma<-array(0,dim=dims)
    dadk <- dat1[seq.int(1, n.row, 3),,,drop=FALSE]
    momk <- dat1[seq.int(2, n.row, 3),,,drop=FALSE]
    kidk <- dat1[seq.int(3, n.row, 3),,,drop=FALSE]

    het <- (mom == 1)
    hethom <- het & (dad == 2)
    ind212<-hethom & (kid == 2)
    ind211<-hethom & (kid == 1)
    # n212 <- colSums(ind212, na.rm=TRUE)
    # n211 <- colSums(ind211, na.rm=TRUE)
    Z.pa[ind212]<-kidk[ind212]-dadk[ind212]
    Z.ma[ind212]<-kidk[ind212]-momk[ind212]
    Z.pa[ind211]<-kidk[ind211]-momk[ind211]
    Z.ma[ind211]<-kidk[ind211]-dadk[ind211]

    hethom <- het & (dad == 0)
    ind011<-hethom & (kid == 1)
    ind010<-hethom & (kid == 0)
    # n011 <- colSums(ind011, na.rm=TRUE)
    # n010 <- colSums(ind010, na.rm=TRUE)
    Z.pa[ind011]<-kidk[ind011]-momk[ind011]
    Z.ma[ind011]<-kidk[ind011]-dadk[ind011]
    Z.pa[ind010]<-kidk[ind010]-dadk[ind010]
    Z.ma[ind010]<-kidk[ind010]-momk[ind010]

    het <- (dad == 1)
    hethom <- het & (mom == 2)
    ind122<-hethom & (kid == 2)
    ind121<-hethom & (kid == 1)
    # n122 <- colSums(ind122, na.rm=TRUE)
    # n121 <- colSums(ind121, na.rm=TRUE)
    Z.pa[ind122]<-kidk[ind122]-dadk[ind122]
    Z.ma[ind122]<-kidk[ind122]-momk[ind122]
    Z.pa[ind121]<-kidk[ind121]-momk[ind121]
    Z.ma[ind121]<-kidk[ind121]-dadk[ind121]

    hethom <- het & (mom == 0)
    ind101<-hethom & (kid == 1)
    ind100<-hethom & (kid == 0)
    # n101 <- colSums(ind101, na.rm=TRUE)
    # n100 <- colSums(ind100, na.rm=TRUE)
    Z.pa[ind101]<-kidk[ind101]-momk[ind101]
    Z.ma[ind101]<-kidk[ind101]-dadk[ind101]
    Z.pa[ind100]<-kidk[ind100]-dadk[ind100]
    Z.ma[ind100]<-kidk[ind100]-momk[ind100]

    het <- (mom == 1) & (dad == 1)
    ind112<-het & (kid == 2)
    ind110<-het & (kid == 0)
    # n112 <- colSums(ind112, na.rm=TRUE)
    # n110 <- colSums(ind110, na.rm=TRUE)
    Z.pa[ind112]<-kidk[ind112]-momk[ind112]
    Z.ma[ind112]<-kidk[ind112]-dadk[ind112]
    Z.pa[ind110]<-kidk[ind110]-dadk[ind110]
    Z.ma[ind110]<-kidk[ind110]-momk[ind110]

    if (xchr){
      ind <- (dad == 1) & (mom == 1) & (kid==1)
      ind0<-sweep(ind,1,sex==0,'&')
      Z.pa[ind0]<-dadk[ind0]
      Z.ma[ind0]<-(-1)*momk[ind0]
      ind0<-sweep(ind,1,sex==1,'&')
      Z.pa[ind0]<-(-1)*dadk[ind0]
      Z.ma[ind0]<-momk[ind0]
    }

    #mat<-cbind(n100,n110,n121,n101,n112,n122,n010,n110,n211,n011,n112,n212)
    #colnames(mat)<-c("p-1","p-1","p-1","p+1","p+1","p+1","m-1","m-1","m-1","m+1","m+1","m+1") #contributions
    out<-list()
    out$Z.pa<-Z.pa
    out$Z.ma<-Z.ma
    return(out)
  }
}
xcontribution <- function(dat,dat1=NULL,dosage=FALSE,sex=NULL){
  #sex: gender info for the offspring, 0 for females and 1 for males
  #Only knockoffs of mothers are required
  n.row <- nrow(dat)
  dad <- dat[seq.int(1, n.row, 3),, drop=FALSE]
  mom <- dat[seq.int(2, n.row, 3),, drop=FALSE]
  kid <- dat[seq.int(3, n.row, 3),, drop=FALSE]
  if (dosage==FALSE){
    Z<-array(0,dim=c(n.row/3,ncol(dat)))
    het <- (mom == 1)
    hethom <- het & (dad == 0)
    Z[hethom & (kid == 0)]<-(-0.5)
    Z[hethom & (kid == 1)]<-0.5
    hethom <- het & (dad == 1)
    ind111<-hethom & (kid == 1)
    Z[hethom & (kid == 0)]<-(-0.5) #110 trios on chrX can only have male offspring
    Z[sweep(ind111,1,sex==0,'&')]<-(-0.5)
    Z[sweep(ind111,1,sex==1,'&')]<-0.5
    Z[hethom & (kid == 2)]<-0.5 #112 trios on chrX can only have female offspring
  }else{
    dims<-dim(dat1)
    dims[1]<-dims[1]/3
    Z<-array(0,dim=dims)
    dadk <- dat1[seq.int(1, n.row, 3),,,drop=FALSE]
    momk <- dat1[seq.int(2, n.row, 3),,,drop=FALSE]
    kidk <- dat1[seq.int(3, n.row, 3),,,drop=FALSE]

    het<-dad==0 & mom==1
    Z[het]<-kidk[het]-0.5*(dadk[het]+momk[het])
    het<-dad==1 & mom==1
    het0<-sweep(het,1,sex==0,'&')
    Z[het0]<-kidk[het0]-0.75*(dadk[het0]+momk[het0])
    het0<-sweep(het,1,sex==1,'&')
    Z[het0]<-kidk[het0]-0.25*(dadk[het0]+momk[het0])
  }
  return(Z)
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
getslidingwindow<-function(left,right,size,pos){
  if (size<=0) stop("The size of sliding window must be a positive number.")
  start<-end<-n<-ind3<-ind4<-c()
  step<-ceiling(size/2) #move half of the window size each time
  l<-left-step
  r<-l+size-1
  while (r<=(right+step)){
    start<-c(start,l)
    end<-c(end,r)
    tmp<-which(pos>=l & pos<=r)
    n<-c(n,length(tmp))
    if (length(tmp)==0) tmp<-0
    ind3<-c(ind3,tmp[1])
    ind4<-c(ind4,tmp[length(tmp)])
    l<-l+step
    r<-r+step
  }
  return(data.frame(start,end,ind3,ind4,n))
}
ACAT<-function(p){
  p[p>0.99]<-0.99 #p[p>1-1e-2]<-1-1e-2
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))

  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}

