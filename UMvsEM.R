library(MASS)
library(pracma)
library(splines)
library(survival)
library(ICsurv)

source("functions.R")

n.variable <- 3
size <- 500
n.spline <- floor(size^(1/3)) + 3
run <- 100
lambda<- 1
spline.ord <- 4

etamatrix <- matrix(0,(n.variable),run)
ICsurvmatrix <- matrix(0,(n.variable),run)

t.seq<-seq(0.1, 4.7, 0.1)
survmatrix.um <- matrix(0, 47, run)
survmatrix.ICsurv <- matrix(0, 47, run)

sumpt1 <- 0
sumpt2 <- 0

for (run_index in 1 : run)	{
  
  cat("\n", "\n", "\n", "run_index=",run_index) 
  
  
  z1<-runif(size, 0, 1)		  
  z2<-rnorm(size, 0, 1)
  z3<-rbinom(size,1,0.5)
  z<-rbind(z1,z2,z3)
  
  beta0<-matrix(c(-1, 0.5, 1.5), n.variable, 1)
  
  u<-runif(size, 0 ,1)
  
  # simulating proportional hazard survival data with exponential baseline survival function
  # \Lamda(t)=t^{1/2} S(t|z)=exp(-t^{1/2}exp(\beta z))
  t<-(log(u) * (-1 / lambda) * c(exp(-t(beta0) %*% z)))^2
  
  delta1 <- rep(0, size)		
  delta2 <- rep(0, size)
  delta3 <- rep(0, size)		  
  
  ctu <- rep(0, size)
  ctv <- rep(0 ,size)
  
  for (ind in 1 : size) {
    
    preuv <- rexp(2, 0.5)
    if (t[ind] <= 5)
    {
      if (t[ind] <= min(preuv)) {
        delta1[ind] <- 1
        ctu[ind] <- min(min(preuv), 5)
        ctv[ind] <- ctu[ind] + 10^(-8)
      } else {
        if (t[ind] > min(preuv) & t[ind] <= max(preuv)) {
          delta2[ind] <- 1
          ctu[ind] <- max(preuv[preuv < t[ind]])
          ctv[ind] <- min(min(preuv[preuv >= t[ind]]), 5)
        } else {										
          if (t[ind] > max(preuv)) {
            delta3[ind] <- 1
            ctv[ind] <- max(preuv)
            ctu[ind] <- ctv[ind] - 10^(-8)
          }
        }			
      }	
    } else {
      delta3[ind] <- 1
      ctv[ind] <- min(max(preuv), 5)
      ctu[ind] <- ctv[ind] - 10^(-8)
    }	
    
  }	
  
  t0 <- proc.time() 
  
  ####################### start the proposed unconstrained maximization method ####################
  knotb <- get_knots(ctu, ctv, delta1, delta2, delta3)
  
  ### bspline results		 
  bspu <- t(splineDesign(knots = knotb, x = ctu, ord = spline.ord))
  ### from bspline to ispline  
  pre.ispu <- apply(bspu, 2, rev)
  
  ispu <- apply(pre.ispu, 2, cumsum)
  ispu <- apply(ispu, 2 ,rev)
  ispu <- ispu[-1,]
  
  ### bspline results
  bspv <- t(splineDesign(knots = knotb, x = ctv, ord = spline.ord))
  ### from bspline to ispline  
  pre.ispv <- apply(bspv, 2, rev)
  
  ispv <- apply(pre.ispv, 2, cumsum)
  ispv <- apply(ispv, 2, rev)
  ispv <- ispv[-1,]
  
  ## create vecters of 1s with length of the numbers of left, interval and right censored###
  number.ones1 <- matrix(1, 1, sum(delta1))
  number.ones2 <- matrix(1, 1, sum(delta2))
  number.ones3 <- matrix(1, 1, sum(delta3))
  
  oldz <- z	  
  index <- 0			
  if (index == 0) {
    z <- oldz
    fixsubx <- matrix(0, 1, size)
    n.total <- n.spline + n.variable		 
  }	 
  
  ############## unconstrained R function for finding minimum point ##########
  out <- optim(par = rep(0, (n.variable + n.spline)), fn = loglike, gr = gradient, method = "BFGS", control = list(reltol = 1e-20))
  
  ### cumulative hazard estimation
  bsp <- t(splineDesign(knots = knotb, x = t.seq, ord = spline.ord, outer.ok = TRUE))
  ### from bspline to ispline  
  pre.isp <- apply(bsp, 2, rev)
  isp <- apply(pre.isp, 2, cumsum)
  isp <- apply(isp, 2 ,rev)
  isp <- isp[-1,]
  
  hazest <- c(matrix(exp(out$par[1 : n.spline]), 1, n.spline) %*% isp)
  survmatrix.um[,run_index] <- exp(-hazest)
  
  etamatrix[,run_index]<-out$par[(n.spline+1) : (n.spline+n.variable)]
  
  #processing time for the proposed unconstrained optimization 
  pt1 <- proc.time() - t0
  sumpt1 <- sumpt1 + pt1
  #############################################################################
  
  ######### estimation using ICsurv package ######################
  t0 <- proc.time()
  
  d1<-rep(0, size)
  d2<-rep(0, size)
  d3<-rep(0, size)
  d1[delta1 == 1] <- 1
  d2[delta2 == 1] <- 1
  d3[delta3 == 1] <- 1

  Li <- rep(0, size)
  Ri <- rep(0, size)
  Li[delta1 == 1] <- 0
  Li[delta2 == 1] <- ctu[delta2 == 1]
  Li[delta3 == 1] <- ctv[delta3 == 1]
  Ri[delta1 == 1] <- ctu[delta1 == 1]
  Ri[delta2 == 1] <- ctv[delta2 == 1]
  Ri[delta3 == 1] <- 0
  
  Xp <- t(z)	
  
  fitsemi <- ICsurv.EM(d1, d2, d3, Li, Ri, Xp, n.int = n.spline-3, order = spline.ord-1, g0 = rep(1,n.spline), b0 = rep(0,n.variable), t.seq = seq(0.1, 4.7, 0.1), tol = 0.001)
  
  survmatrix.ICsurv[,run_index] <- exp(-c(fitsemi$hz))
  ICsurvmatrix[,run_index]<-fitsemi$b		
  
  # processing time for the ICsurv EM algorithm 
  pt2 <- proc.time() - t0
  sumpt2 <- sumpt2 + pt2
  #####################################################
}

etamean <- rowMeans(etamatrix)
c(beta0)
etamean
etamean - c(beta0)
ICmean <- rowMeans(ICsurvmatrix)
ICmean
ICmean - c(beta0)

sumpt1
sumpt2
sumpt2/sumpt1

surv.true <- exp(-t.seq^(1/2))
meansurv.um <- rowMeans(survmatrix.um)
meansurv.ICsurv <- rowMeans(survmatrix.ICsurv)

pdf("survivals.pdf")	
plot(c(-0.2, 5), c(-0.05, 1.05),type = "n", xlab = " ", ylab = " ")
lines(t.seq, surv.true, lty = 1, col = 1, lwd = 2)
lines(t.seq, meansurv.um, lty = 2, col = 2, lwd = 2)
lines(t.seq, meansurv.ICsurv, lty = 4, col = 4, lwd = 2)
legend(2.5, 0.8, c("True", "UM", "EM"), lty = c(1,2,4), col = c(1,2,4))
dev.off()
