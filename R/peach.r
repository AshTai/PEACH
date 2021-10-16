#' Detection of cell separation–induced gene expression through a penalized deconvolution approach.
#'
#' PEACH adopts a deconvolution method to detect cell separation–induced genes.
#'
#' @name PEACH
#' @author An-Shun Tai \email{anshuntai@nctu.edu.tw}
#' @param res.set A G by N expression matrix of reference, where G is the gene number and N is the total number of reference.
#' @param mix.set A G by M expression matrix of bulk samples, where G is the gene number and M is the total number of bulk samples.
#' @param iter The number of iterations. Default is 10000.
#' @param parC A parameter of the spike-and-slab distribution. Default is 10^-4.
#' @return A list of estimated cell proportions (w),
#' the probability of informative genes (p),
#' the fraction of informative genes per cell type (pi), and
#' the cell-specific profiles (mean).
#' @export
#' @examples


peach_main <- function(ref,bulk,iter=10^4,parC=10^-4){
  require(MCMCpack)
  require(invgamma)

  n.type <- ncol(ref)
  n.gene <- nrow(ref)

  hyper.par <- list(eta2=10^6,c=parC,c0=1,d0=1)
  prior.t <- prior.gamma <- matrix(NA,iter,n.gene)
  prior.pi <- matrix(NA,iter,1)
  prior.w <- matrix(NA,iter,n.type)

  nL.w <- function(x,y,z){
    w.m <- x
    value.nL <- mean((y-z%*%matrix(w.m))^2)
    value.nL
  }
  optim.w <- function(y,z){
    ff <- constrOptim(rep(1/(n.type),n.type), f=nL.w,grad=NULL,y=y,z=z,ui = rbind(rep(1,n.type),-rep(1,n.type), diag(n.type)), ci = c(0.99,-1.01,rep(0,n.type)))
    opt.w <- ff$par
    c(opt.w)
  }
  lm.w <- optim.w(y=bulk,z=ref)
  lm.w <- lm.w/sum(lm.w)

  lm.res <- bulk - ref%*%matrix(lm.w)
  t.ini.ord <- order(abs(lm.res),decreasing = T)[1:300]
  t.ini <- rep(0,n.gene); t.ini[t.ini.ord] <- 1
  gamma.ini <- lm.res

  prior <- list(
    gamma=gamma.ini,
    w=lm.w,
    t=t.ini,
    pi=0.1,
    sigma2=mean(lm.res^2)
  )

  iter.c <- 1; stop.sign <- F; count <- 1
  run.time <- Sys.time()
  while(stop.sign==F){
    if(count==floor(iter.c/iter*100)){
      tt <- difftime(Sys.time(),run.time,units ="mins")
      print(paste(count,"%","--","Expected running time =",round((100-count)*tt,3),"mins"))  ;
      count <- count +1
      run.time <- Sys.time()
    }

    gamma.sampling <- function(d){
      w1 <- 1/prior$sigma2; w2 <- 1/hyper.par$eta2*prior$t+1/(hyper.par$eta2*hyper.par$c)*(1-prior$t)
      pos.mean <- w1/(w1+w2)*(bulk-ref%*%matrix(prior$w))
      pos.var <- 1/(w1+w2)
      rnorm(n.gene,pos.mean,sqrt(pos.var))
    }
    prior$gamma <- gamma.sampling(1)
    prior.gamma[iter.c,] <- prior$gamma

    # update sigma2
    rss <- sum((bulk-ref%*%matrix(prior$w)-prior$gamma)^2)
    prior$sigma2 <- rinvgamma(1,(n.gene)/2-1,(rss)/2)

    # sampling t
    density1 <- dnorm(prior$gamma,0,sqrt(hyper.par$eta2),log = T) +
      log(prior$pi)
    density2 <-  dnorm(prior$gamma,0,sqrt(hyper.par$eta2*hyper.par$c),log = T) +
      log((1-prior$pi))
    prob.gamma <-  1/(1+exp(density2-density1))
    prior$t <-  rbinom(n=n.gene,1,prob=prob.gamma)
    prior.t[iter.c,] <- prior$t

    # update pi
    prior$pi <- rbeta(1,hyper.par$c0+sum(prior$t==1),
                      hyper.par$d0+n.gene-sum(prior$t==1))
    prior.pi[iter.c] <- prior$pi


    # update component prop
    gene_index <- which(prior$t==0)
    w_s <- rdirichlet(1,alpha =prior$w*100+1)

    mu.mix.1 <- ref%*%matrix(w_s)
    mu.mix.2 <- ref%*%matrix(prior$w)

    density1 <- (dnorm(bulk,mean=mu.mix.1,sd=sqrt(prior$sigma2),log = T)[gene_index])
    density2 <- (dnorm(bulk,mean=mu.mix.2,sd=sqrt(prior$sigma2),log = T)[gene_index])
    density1[density1==-Inf] <- 0; density2[density2==-Inf] <- 0
    r1 <- sum(density1)
    r2 <- sum(density2)
    if(exp(r1-r2)>runif(1)){
      prior$w <- w_s
    }
    prior.w[iter.c,] <- prior$w

    iter.c <- iter.c+1
    if(iter.c > iter){stop.sign <- T}
  }#end while
  structure(list(w=prior.w,t=prior.t,pi=prior.pi,gamma=prior.gamma))
}


peach <- function(ref,bulk,iter=10^4,seqC=10^-c(1:6)){
  res.w <- matrix(NA,ncol(ref),length(seqC))
  res.t <- matrix(NA,nrow(ref),length(seqC))
  res.gamma <- matrix(NA,nrow(ref),length(seqC))

  for(k in 1:length(seqC)){
    f <- peach_main(ref=ref,bulk=bulk,iter=iter,parC=seqC[k])
    res.w[,k] <- colMeans(f$w[round(iter*0.7):iter,])
    res.t[,k] <- colMeans(f$t[round(iter*0.7):iter,])
    res.gamma[,k] <- colMeans(f$gamma[round(iter*0.7):iter,])
    res.gamma[,k][res.t[,k]<0.5] <- 0
  }
  loss <- colSums((bulk%*%matrix(1,1,length(seqC))-ref%*%res.w-res.gamma)^2)
  structure(list(w=res.w,t=res.t,gamma=res.gamma,loss=loss))
}






