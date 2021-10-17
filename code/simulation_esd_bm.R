################################################################################
## Simulation - spectral density of the block sthocastic model
################################################################################


# Setting working directory ----------------------------------------------------

wd <- "."

setwd(wd)


# Dependecies ------------------------------------------------------------------

library('statGraph')
library('igraph')
library(snowfall)
source('statGraph.R')


# Parameters -------------------------------------------------------------------

executionName <- "esd_BM"

scenarios <- 1:5

# Generation of BM graphs ------------------------------------------------------

BM <- function(sizes, p0, p) {
  M <- length(sizes)
  n <- sum(sizes)
  theta <- matrix(0,M,M)
  #s <- sum(upper.tri(theta))
  theta[upper.tri(theta)] <- p0
  theta <- theta + t(theta)
  diag(theta) <- p
  g <-sample_sbm(n, pref.matrix = theta, block.sizes = sizes)
  return(as_adjacency_matrix(g))
}


# BM theoretical density ------------------------------------------------------

# Eq(22) from Avrachenkov et al corresponds to c.z=F(c.z,z). Here we code function F: Complex \to Complex
# The function can be applied on vectors of values as follows.
# arguments c.z is a list with $re and $im. Both are matrices with M columns and as many rows as points z. 
#           (x,y) are real and imaginary part of points z. It is a list of vectors with as many entries as points z.
#           NB x,y should be of same length. as well as the number of rows or real and imaginary part of c.z
# Returns a list with $re and $im. Both are matrices with M columns and as many rows as points z.
EQ.22.old <- function(c.z,x,y,zeta0,zeta){
  M <- length(zeta)
  if (length(x)==1){
    real <- rep(NA,M)
    imaginary <- rep(NA,M)
    for (m in 1:M){
      sq.module <- (x+zeta[m]*c.z$re[m] +zeta0*(sum(c.z$re)-c.z$re[m]))^2 + (y+zeta[m]*c.z$im[m] +zeta0*(sum(c.z$im)-c.z$im[m]))^2
      real[m] <- -(x+zeta[m]*c.z$re[m] +zeta0*(sum(c.z$re)-c.z$re[m]))/sq.module/M
      imaginary[m] <- (y+zeta[m]*c.z$im[m] +zeta0*(sum(c.z$im)-c.z$im[m])) /sq.module/M
    }
    
  }
  else {
    real <- matrix(NA,length(x),M)
    imaginary <- matrix(NA,length(y),M)
    for (m in 1:M){
      sq.module <- (x+zeta[m]*c.z$re[,m] +zeta0*(rowSums(c.z$re)-c.z$re[,m]))^2 + (y+zeta[m]*c.z$im[,m] +zeta0*(rowSums(c.z$im)-c.z$im[,m]))^2
      real[,m] <- -(x+zeta[m]*c.z$re[,m] +zeta0*(rowSums(c.z$re)-c.z$re[,m]))/sq.module/M
      imaginary[,m] <- (y+zeta[m]*c.z$im[,m] +zeta0*(rowSums(c.z$im)-c.z$im[,m])) /sq.module/M
    }
  }
  return(list(re=real, im=imaginary))
}

EQ.22 <- function(c.z,x,y,zeta0,zeta){
  M <- length(zeta)
  if (length(zeta0) == 1) {
    zeta0 <- matrix(zeta0, M, M)
    diag(zeta0) <- 0
  }
  else {
    tmp <- matrix(0, M, M)
    tmp[upper.diag(tmp)] <- zeta0
    zeta0 <- tmp + t(tmp)
  }
  if (length(x)==1){
    real <- rep(NA,M)
    imaginary <- rep(NA,M)
    for (m in 1:M){
      sq.module <- (x+zeta[m]*c.z$re[m] + sum(zeta0[m,]*(c.z$re)))^2 + (y+zeta[m]*c.z$im[m] + sum(zeta0[m,]*(c.z$im)))^2
      real[m] <- -(x+zeta[m]*c.z$re[m] + sum(zeta0[m,]*(c.z$re)))/sq.module/M
      imaginary[m] <- (y+zeta[m]*c.z$im[m] + sum(zeta0[m,]*(c.z$im))) /sq.module/M
    }
    
  }
  else {
    real <- matrix(NA,length(x),M)
    imaginary <- matrix(NA,length(y),M)
    for (m in 1:M){
      sq.module <- (x+zeta[m]*c.z$re[,m] + rowSums(apply(c.z$re,2,'*',zeta0[m,])))^2 + (y+zeta[m]*c.z$im[,m] +rowSums(apply(c.z$im,2,'*',zeta0[m,])))^2
      real[,m] <- -(x+zeta[m]*c.z$re[,m] + rowSums(apply(c.z$re,2,'*',zeta0[m,])))/sq.module/M
      imaginary[,m] <- (y+zeta[m]*c.z$im[,m] + rowSums(apply(c.z$im,2,'*',zeta0[m,]))) /sq.module/M
    }
  }
  return(list(re=real, im=imaginary))
}

# x=c(1,2)
# y=c(2,3)
# c.z=list(re=matrix(rep(x,M),nrow=length(x)),im=matrix(rep(y,M),nrow=length(y)))
# EQ.22(c.z,x,y,zeta)
# x=2
# y=3
# c.z=list(re=rep(x,M),im=rep(y,M))
# EQ.22(c.z,x,y,zeta)

# Compute the Stieltjes transform of the limit. Corollary 3.
# Argument is a complex number z=x+iy. 
# Returns a complex number in the form of a real part and an imaginary part 
stieltjes <- function(x,y,zeta0,zeta,eps=1e-6){
  M <- length(zeta)
  # How to choose good initial value for c.z ???? I use z...
  c.z <- list(re=matrix(rep(x,M),nrow=length(x)),im=matrix(rep(y,M),nrow=length(y)))  
  not.converged <- TRUE
  # search for the solution of Eq (22)
  while (not.converged){
    c.z.old <- c.z
    c.z <- EQ.22(c.z,x,y,zeta0,zeta) 
    # Here maybe I should include the condition (23) on imaginary parts...
    not.converged <- max(c(abs(c.z.old$re -c.z$re),abs(c.z.old$im -c.z$im))) > eps
  }
  res <- if (length(x)==1) {list(re=sum(c.z$re),im=sum(c.z$im))} else list(re=apply(c.z$re,1,sum),im=apply(c.z$im,1,sum))
  return(res)
}

# Inverse of the stieljes function. Recover the density from its Stieltjes transform 
# I found this formula on the wiki webpage
# https://en.wikipedia.org/wiki/Stieltjes_transformation
# Beware: Avrachenkov's definition and wiki's def of Stieltjes transform differ by a minus sign
inv.stieltjes <- function(x,zeta0,zeta,eps=1e-6){
  return((stieltjes(x,rep(eps,length(x)),zeta0,zeta)$im - stieltjes(x,rep(-eps,length(x)),zeta0,zeta)$im)/2/pi)
}


# Centered BM adjancency matrix -----------------------------------------------

centered.BM <- function(A, sizes, p0, p) {
  Amean <- matrix(p0, nrow(A), ncol(A))
  M <- length(sizes)
  for (i in 1:M) {
    from <- 1
    if (i > 1)
      from <- sum(sizes[1:(i-1)])
    to <- sum(sizes[1:i])
    Amean[from:to, from:to] <- p[i]
  }

  return(Amean)
}


# Plotting densities -----------------------------------------------------------

plot.BM.densities <- function(esd, fs, sizes, p0, p, cols=1:length(fs)) {
  theoretical <- list()
  for (i in 1:length(fs)) {
    theoretical[[i]] <- fs[[i]](esd$x, sizes, p0, p)
  }
  plot.densities(esd, theoretical, cols)
}

plot.densities <- function(esd, theoretical, cols=1:length(theoretical)) {
  f <- list()
  Min <- Inf
  Max <- -Inf
  for (i in 1:length(theoretical)) {
    Min <- min(Min, theoretical[[i]])
    Max <- max(Max, theoretical[[i]])
  }
  Min <- min(esd$y, Min)
  Max <- max(esd$y, Max)
  plot(esd, ylim=c(Min, Max))
  for (i in 1:length(theoretical)) {
    lines(esd$x, theoretical[[i]],col=cols[i], lwd=3)
  }
}

# Scenarios --------------------------------------------------------------------

run <- function(scenario) {
  setwd(wd)

  system(paste("renice 10 -p", Sys.getpid()))

  block.sizes <- c(100, 80, 300, 150, 50, 200, 75, 180, 70, 90)
  p <-c(0.8, 0.5, 0.6, 0.7, 0.4, 0.9, 0.55, 0.42, 0.38, 0.86) # intra-group proba
  p0 <- 0.2
  M <- 3
  K <- 300

  # Avra scenario, with 3 communities
  if (scenario == 1) {
      block.sizes <- rep(K,M)
  }
  # Avra scenario, with more communities
  if (scenario == 2) {
      M <- 10
      block.sizes <- rep(K,M)
  }
  # Avra scenario with different block sizes 
  if (scenario == 3) {
      block.sizes <- block.sizes[1:M]  
  }
  # Avra scenario with larger p0
  if (scenario == 4) {
      block.sizes <- rep(K,M) 
      p0 <- 0.9
  }
  # Avra scenario with different p0
  if (scenario == 5) {
      p0 <- c(0.1, 0.2, 0.05)
      block.sizes <- rep(K,M) 
  }
  # Avra scenario with different p0
  if (scenario == 6) {
      M <- 5
      p0 <- c(0.1, 0.2, 0.7, 0.3, 0.21, 0.36, 0.08, 0.09, 0.25, 0.5)
      block.sizes <- rep(K,M) 
  }
  n <- sum(block.sizes)
  p <- p[1:M]
  pmax <- max(p)
  v <- n*pmax*(1-pmax)
  #str(J.K)
  A <- BM(block.sizes, p0, p)
  #pmax2 <- max(pmax, p0)
  #v2 <- ncol(A)*pmax2*(1-pmax2)
  A.mean <- centered.BM(A, block.sizes, p0, p) # here the groups are ordered
  #str(A.mean)
  A.norm <- (A-A.mean)/sqrt(v)
  esd.norm <- spectralDensity(A.norm*sqrt(n)) # A.norm is already scaled

  # Limit of normalized ESD - according to Avrachenkov
  zeta <- p*(1-p)/pmax/(1-pmax)
  zeta0 <- p0*(1-p0)/pmax/(1-pmax)

  # Limit of normalized ESD - according to Avrachenkov
  #zeta2 <- p*(1-p)/pmax2/(1-pmax2)
  #zeta02 <- p0*(1-p0)/pmax2/(1-pmax2)

  #str(A.mean)

  A.norm <- (A-A.mean)/sqrt(v)
  esd.norm <- spectralDensity(A.norm*sqrt(n)) # A.norm is already scaled

  #A.norm2 <- (A-A.mean)/sqrt(v2)
  #esd.norm2 <- spectralDensity(A.norm2*sqrt(n)) # A.norm is already scaled

  esd <- spectralDensity((A/sqrt(v))*sqrt(n))

  x <- esd.norm$x[seq(1,1024,by=4)]
  y <- rep(NA, length(x))
  #y <- inv.stieltjes(x, zeta0, zeta)
  y[1:120] <- inv.stieltjes(x[1:120], zeta0, zeta)
  y[121:130] <- inv.stieltjes(x[121:130], zeta0, zeta)
  result <- list("block.sizes"=block.sizes, "p"=p, "p0"=p0, "K"=K, 
                 "n"=sum(block.sizes), "esd"=esd, "esd.norm"=esd.norm, 
                 "x"=x, "y"=y)
  save(result, file=paste(executionName, "_", scenario, ".RData", sep=""))
  y[135:256] <-inv.stieltjes(x[135:256], zeta0, zeta)
  #difficult to compute for x close to 0
  y[131] <- inv.stieltjes(x[131], zeta0, zeta)
  result <- list("block.sizes"=block.sizes, "p"=p, "p0"=p0, "K"=K, 
                 "n"=sum(block.sizes), "esd"=esd, "esd.norm"=esd.norm, 
                 "x"=x, "y"=y)
  save(result, file=paste(executionName, "_", scenario, ".RData", sep=""))
  y[134] <- inv.stieltjes(x[134], zeta0, zeta)
  y[133] <- inv.stieltjes(x[133], zeta0, zeta)
  y[132] <- inv.stieltjes(x[132], zeta0, zeta)

  result <- list("block.sizes"=block.sizes, "p"=p, "p0"=p0, "K"=K, 
                 "n"=sum(block.sizes), "esd"=esd, "esd.norm"=esd.norm, 
                 "x"=x, "y"=y)
  save(result, file=paste(executionName, "_", scenario, ".RData", sep=""))
  return(result)
} 

# Setting up parallel execution ------------------------------------------------

sfInit(parallel = TRUE , type ="SOCK", socketHosts =
       c(rep("localhost", length(scenarios))))

# Exporting global variables ---------------------------------------------------

sfExportAll() # export all variables
sfLibrary(igraph) # export packages


# Simulation experiments -------------------------------------------------------

results <- sfClusterApplySR(scenarios, run, name=executionName, restore=F, perUpdate=1)

sfStop()

save(results, file=paste(executionName, ".RData", sep=""))

