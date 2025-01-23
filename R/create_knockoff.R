#' Create knockoff genotype data
#'
#' Create knockoff genotype data using phased haplotype data.
#'
#' @param trio.hap A 6n*p matrix for trio haplotype data, in which n is the number of trios and p is the number of variants. Each trio must consist of father, mother, and offspring (in this order). The haplotypes must be coded as 0 or 1. Missing haplotypes are not allowed.
#' @param duo.hap A 4m*p matrix for duo haplotype data, in which m is the number of duos and p is the number of variants. Each duo must consist of a single parent and offspring (in this order). The haplotypes must be coded as 0 or 1. Missing haplotypes are not allowed.
#' @param pos A numeric vector of length p for the position of p variants.
#' @param M A positive integer for the number of knockoffs. The default is 10.
#' @param maxcor A real number in a range of [0,1] for hierarchical clustering of neighboring variants used to generate knockoff parents. The default is 0.7.
#' @param maxbp A positive integer for the size of neighboring base pairs used to generate knockoff parents. The default is 80000.
#' @param phasing.dad A numeric vector of length n that contains 1 or 2 to indicate which paternal haplotype was transmitted to offspring in each trio. If NULL, the function will calculate the phasing information based on the input trio haplotype matrix.
#' @param phasing.mom A numeric vector of length n that contains 1 or 2 to indicate which maternal haplotype was transmitted to offspring in each trio. If NULL, the function will calculate the phasing information based on the input trio haplotype matrix.
#' @param seed An integer for the random seed used for knockoff generation.
#' @return A list that contains:
#' \describe{
#'   \item{trio.ko}{A 3n*p*M array for knockoff trio genotype data if trio.hap is provided.}
#'   \item{duo.ko}{A 3m*p*M array for knockoff duo genotype data if duo.hap is provided.}
#'   \item{duo}{A 3m*p matrix for duo genotype data if duo.hap is provided.}
#' }
#' @importFrom stats as.dist cutree hclust median pcauchy pnorm princomp
#' @export
#' @examples
#' data(KnockoffTrio.example)
#' knockoff<-create_knockoff(trio.hap=KnockoffTrio.example$trio.hap,
#'           duo.hap=KnockoffTrio.example$duo.hap, pos=KnockoffTrio.example$pos, M=10)
create_knockoff<-function(trio.hap=NULL, duo.hap=NULL, pos, M=10, maxcor=0.7, maxbp=80000, phasing.dad=NULL, phasing.mom=NULL, seed=100){
  set.seed(seed)
  if (is.null(trio.hap) & is.null(duo.hap)) stop("The input trio and duo haplotype matrices cannot both be NULL at the same time.")
  trio.ko<-duo.ko<-duo<-NULL
  if (!is.null(trio.hap)){
    if (nrow(trio.hap) %% 6!=0) stop("The number of rows of the input trio haplotype matrix must be a multiple of six.")
    if (!all(trio.hap %in% c(0,1))) stop("The input haplotype matrix can only contain 0 or 1.")
    n<-nrow(trio.hap)/6 #number of trios
    nsnp<-ncol(trio.hap) #number of variants

    ind.hap.off<-c((1:n)*6-1,(1:n)*6)
    if (is.null(phasing.dad) | is.null(phasing.mom)){
      phasing<-get_phasing(trio.hap)
      phasing.dad<-phasing$phasing.dad
      phasing.mom<-phasing$phasing.mom
    }
    phasing.dad<-phasing.dad%%2
    phasing.mom<-phasing.mom%%2
    #remove offspring
    trio.hap<-trio.hap[-ind.hap.off,,drop=FALSE]

    dad_hap<-trio.hap[sort(c(seq(1,(4*n),4),seq(2,(4*n),4))),,drop=FALSE]
    ind_dad<-2*(1:n)+sample(c(0,-1),n,replace = T)
    ind_dad2<-(1:(2*n))[-ind_dad]

    mom_hap<-trio.hap[sort(c(seq(3,(4*n),4),seq(4,(4*n),4))),,drop=FALSE]
    ind_mom<-2*(1:n)+sample(c(0,-1),n,replace = T)
    ind_mom2<-(1:(2*n))[-ind_mom]

    dadk1<-create_knockoff_parent(dad_hap[ind_dad,,drop=FALSE],pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    dadk2<-create_knockoff_parent(dad_hap[ind_dad2,,drop=FALSE],pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    momk1<-create_knockoff_parent(mom_hap[ind_mom,,drop=FALSE],pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    momk2<-create_knockoff_parent(mom_hap[ind_mom2,,drop=FALSE],pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    dadk<-dadk1+dadk2
    momk<-momk1+momk2
    kidk<-array(dim=dim(dadk))

    for (i in 1:n){
      if (ind_dad[i]%%2==phasing.dad[i]) hap1<-dadk1[i,,,drop=FALSE] else hap1<-dadk2[i,,,drop=FALSE]
      if (ind_mom[i]%%2==phasing.mom[i]) hap2<-momk1[i,,,drop=FALSE] else hap2<-momk2[i,,,drop=FALSE]
      kidk[i,,]<-hap1+hap2
    }

    trio.ko<-array(dim=c(3*n,nsnp,M)) #knockoff for trio, dim=(3*n,nsnp,M)
    for (i in 1:n){
      trio.ko[i*3-2,,]<-dadk[i,,]
      trio.ko[i*3-1,,]<-momk[i,,]
      trio.ko[i*3,,]<-kidk[i,,]
    }
  }
  if (!is.null(duo.hap)){ #data has duos
    if (nrow(duo.hap) %% 4!=0) stop("The number of rows of the input duo haplotype matrix must be a multiple of four.")
    if (!all(duo.hap %in% c(0,1))) stop("The input haplotype matrix can only contain 0 or 1.")
    n<-nrow(duo.hap)/4 #number of duos
    nsnp<-ncol(duo.hap) #number of variants

    ind.hap.off<-sort(c((1:n)*4-1,(1:n)*4))
    phasing<-get_phasing(trio.hap=NULL, duo.hap = duo.hap, xchr = FALSE)
    phasing.duo.parent<-phasing$phasing.duo.parent
    phasing.duo.offspring<-phasing$phasing.duo.offspring

    phasing.duo.parent<-phasing.duo.parent%%2
    phasing.duo.offspring<-phasing.duo.offspring%%2

    #remove offspring
    dad_hap<-duo.hap[-ind.hap.off,,drop=FALSE]
    ind_dad<-2*(1:n)+sample(c(0,-1),n,replace = T)
    ind_dad2<-(1:(2*n))[-ind_dad]

    off_hap<-duo.hap[ind.hap.off,,drop=FALSE]
    ind_mom<-2*(1:n)-phasing.duo.offspring #index for the haplotype transmitted to offspring from the missing parent

    dadk1<-create_knockoff_parent(dad_hap[ind_dad,,drop=FALSE],pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    dadk2<-create_knockoff_parent(dad_hap[ind_dad2,,drop=FALSE],pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    momk1<-create_knockoff_parent(off_hap[ind_mom,,drop=FALSE],pos,M,corr_max=maxcor,maxBP.neighbor=maxbp)
    if (is.null(trio.hap)) maf<-colMeans(dad_hap, na.rm = TRUE) else{
      maf1<-colMeans(trio.hap, na.rm = TRUE) #MAF for parents in trios; no need to "/2" because it is haplotype
      maf2<-colMeans(dad_hap, na.rm = TRUE) #MAF for single parents in duos; no need to "/2" because it is haplotype
      n1<-nrow(trio.hap)
      n2<-nrow(dad_hap)
      maf<-(maf1*n1+maf2*n2)/(n1+n2)
    }
    momk2<-array(dim=dim(momk1))
    for (i in 1:dim(momk2)[1]) momk2[i,,]<-maf #E(knockoff missing haplotype)=MAF

    dadk<-dadk1+dadk2
    momk<-momk1+momk2
    kidk<-array(dim=dim(dadk))

    for (i in 1:n){
      if (ind_dad[i]%%2==phasing.duo.parent[i]) hap1<-dadk1[i,,,drop=FALSE] else hap1<-dadk2[i,,,drop=FALSE]
      kidk[i,,]<-hap1+momk1[i,,,drop=FALSE]
    }

    duo.ko<-array(dim=c(3*n,nsnp,M)) #knockoff genotype for duos in the form of trios, dim=(3*n,nsnp,M); missing parent = transmitted haplotype to offspring + MAF
    for (i in 1:n){
      duo.ko[i*3-2,,]<-dadk[i,,]
      duo.ko[i*3-1,,]<-momk[i,,]
      duo.ko[i*3,,]<-kidk[i,,]
    }

    duo<-array(dim=c(3*n,nsnp)) #original genotype for duos in the form of trios, dim=(3*n,nsnp); missing parent = transmitted haplotype to offspring + MAF
    for (i in 1:n){
      duo[i*3-2,]<-dad_hap[2*i-1,]+dad_hap[2*i,]
      duo[i*3-1,]<-off_hap[2*i-phasing.duo.offspring[i],]+maf
      duo[i*3,]<-off_hap[2*i-1,]+off_hap[2*i,]
    }
  }
  out<-list()
  out$trio.ko<-trio.ko #knockoff genotype for trios
  out$duo.ko<-duo.ko #knockoff genotype for duos (in the form of trios)
  out$duo<-duo #original genotype for duos (in the form of trios)
  return(out)
}
create_knockoff_parent<-function(X,pos,M=10,corr_max=0.75,maxN.neighbor=Inf,maxBP.neighbor=100000) {

  X <- as.matrix(X)
  sparse.fit <- sparse.cor(X)
  cor.X <- sparse.fit$cor
  cor.X[is.na(cor.X)] <- 0
  cor.X[is.infinite(cor.X)] <- 0
  cov.X <- sparse.fit$cov
  cov.X[is.na(cov.X)] <- 0
  cov.X[is.infinite(cov.X)] <- 0
  Sigma.distance = as.dist(1 - abs(cor.X))
  if(ncol(X)>1){
    fit = hclust(Sigma.distance, method="single")
    corr_max = corr_max
    clusters = cutree(fit, h=1-corr_max)
  }else{clusters<-1}

  X_k<-array(0,dim=c(nrow(X),ncol(X),M));index.exist<-c()
  for (k in unique(clusters)){
    cluster.fitted<-cluster.residuals<-matrix(NA,nrow(X),sum(clusters==k))
    indk<-which(clusters==k)
    for(i in indk){
      index.pos<-which(pos>=max(pos[i]-maxBP.neighbor,pos[1]) & pos<=min(pos[i]+maxBP.neighbor,pos[length(pos)]))
      temp<-abs(cor.X[i,]);temp[indk]<-0;temp[-index.pos]<-0

      index<-order(temp,decreasing=T)
      index<-setdiff(index[1:min(length(index),sum(temp>0.05),floor((nrow(X))^(1/3)),maxN.neighbor)],i)

      y<-X[,i]
      if(length(index)==0){fitted.values<-mean(y)}else{

        x<-X[,index,drop=F];temp.xy<-rbind(mean(y),crossprod(x,y)/length(y)-colMeans(x)*mean(y))
        x.exist<-c()
        for(j in 1:M){
          x.exist<-cbind(x.exist,X_k[,intersect(index,index.exist),j])
        }
        temp.xy<-rbind(temp.xy,crossprod(x.exist,y)/length(y)-colMeans(x.exist)*mean(y))

        temp.cov.cross<-sparse.cov.cross(x,x.exist)$cov
        temp.cov<-sparse.cor(x.exist)$cov
        temp.xx<-cov.X[index,index]
        temp.xx<-rbind(cbind(temp.xx,temp.cov.cross),cbind(t(temp.cov.cross),temp.cov))

        temp.xx<-cbind(0,temp.xx)
        temp.xx<-rbind(c(1,rep(0,ncol(temp.xx)-1)),temp.xx)

        pca.fit<-princomp(covmat=temp.xx)
        v<-pca.fit$loadings
        cump<-cumsum(pca.fit$sdev^2)/sum(pca.fit$sdev^2)
        n.pc<-which(cump>=0.999)[1]
        pca.index<-intersect(1:n.pc,which(pca.fit$sdev!=0))

        temp.inv<-v[,pca.index,drop=F]%*%(pca.fit$sdev[pca.index]^(-2)*t(v[,pca.index,drop=F]))
        temp.beta<-temp.inv%*%temp.xy

        temp.j<-1
        fitted.values<-temp.beta[1]+crossprod(t(x),temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])-sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),,drop=F])
        temp.j<-temp.j+ncol(x)
        for(j in 1:M){
          temp.x<-as.matrix(X_k[,intersect(index,index.exist),j])
          if(ncol(temp.x)>=1){
            fitted.values<-fitted.values+crossprod(t(temp.x),temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])-sum(colMeans(temp.x)*temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),,drop=F])
          }
          temp.j<-temp.j+ncol(temp.x)
        }
      }
      residuals<-y-fitted.values
      cluster.fitted[,match(i,indk)]<-as.vector(fitted.values)
      cluster.residuals[,match(i,indk)]<-as.vector(residuals)

      index.exist<-c(index.exist,i)
    }
    cluster.sample.index<-sapply(1:M,function(x)sample(1:nrow(X)))
    for(j in 1:M){
      X_k[,indk,j]<-cluster.fitted+cluster.residuals[cluster.sample.index[,j],,drop=F]
    }
  }
  return(X_k)
}
get_phasing<-function(trio.hap, duo.hap=NULL, xchr=FALSE){
  #Per SHAPEIT's duoHMM, a child's haplotypes are now ordered with the paternal haplotype first and the maternal haplotype second.
  #for X chr, female offsprings haplotypes may be ordered with maternal first and paternal second
  out<-list()
  if (!is.null(trio.hap)){
    nsnp<-ncol(trio.hap)
    ntrio<-nrow(trio.hap)/6
    phasing.dad<-phasing.mom<-rep(1,ntrio)
    if (!xchr){
      for (i in 1:ntrio){
        index<-((i-1)*6+1):(i*6)
        if (length(which(trio.hap[index[1],]==trio.hap[index[5],]))/nsnp>=length(which(trio.hap[index[2],]==trio.hap[index[5],]))/nsnp) phasing.dad[i]<-1 else phasing.dad[i]<-2
        if (length(which(trio.hap[index[3],]==trio.hap[index[6],]))/nsnp>=length(which(trio.hap[index[4],]==trio.hap[index[6],]))/nsnp) phasing.mom[i]<-1 else phasing.mom[i]<-2
      }
    } else{
      #father's phasing info is not important for X chr and depends on proband's sex
      for (i in 1:ntrio){
        index<-((i-1)*6+1):(i*6)
        #first, determine which of the female proband's haptypes was inherited from the father
        momhap<-6
        if (length(which(trio.hap[index[1],]==trio.hap[index[5],]))/nsnp<length(which(trio.hap[index[1],]==trio.hap[index[6],]))/nsnp) momhap<-5
        if (length(which(trio.hap[index[3],]==trio.hap[index[momhap],]))/nsnp>=length(which(trio.hap[index[4],]==trio.hap[index[momhap],]))/nsnp) phasing.mom[i]<-1 else phasing.mom[i]<-2
      }
    }
    out$phasing.dad<-phasing.dad
    out$phasing.mom<-phasing.mom
  }
  if (!is.null(duo.hap)){
    nsnp<-ncol(duo.hap)
    nduo<-nrow(duo.hap)/4
    phasing.duo.parent<-phasing.duo.offspring<-rep(1,nduo)
    for (i in 1:nduo){
      index<-((i-1)*4+1):(i*4)
      percent<-c(length(which(duo.hap[index[1],]==duo.hap[index[3],])),length(which(duo.hap[index[1],]==duo.hap[index[4],])),length(which(duo.hap[index[2],]==duo.hap[index[3],])),length(which(duo.hap[index[2],]==duo.hap[index[4],])))
      ind<-which.max(percent)
      if (ind==1 | ind==2) phasing.duo.parent[i]<-1 else phasing.duo.parent[i]<-2 #which parental haplotype was transmitted from the single parent to offspring
      if (ind==1 | ind==3) phasing.duo.offspring[i]<-2 else phasing.duo.offspring[i]<-1 #which offspring haplotype is from the missing parent
    }
  out$phasing.duo.parent<-phasing.duo.parent
  out$phasing.duo.offspring<-phasing.duo.offspring
  }
  return(out)
}
sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}
sparse.cov.cross <- function(x,y){
  n <- nrow(x)
  cMeans.x <- colMeans(x);cMeans.y <- colMeans(y)
  covmat <- (as.matrix(crossprod(x,y)) - n*tcrossprod(cMeans.x,cMeans.y))/(n-1)
  list(cov=covmat)
}
