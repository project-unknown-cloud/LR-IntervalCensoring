############################ the log likelihood function ######################
loglike <- function(eta)
{
	eta <- matrix(c(eta), (n.total),1)
	alpha<-eta[(1):(n.spline),]
	subu<-c(t(exp(alpha))%*%ispu)
	subv<-c(t(exp(alpha))%*%ispv)
	
	if (n.total > n.spline)	{
		beta<-eta[(n.spline+1):(n.total),]
		subx<-exp(c(t(beta) %*% z + fixsubx))	
	} else
		subx <- exp(c(fixsubx))
		
	log1 <- sum(log(1 - exp(-subu[delta1==1] * subx[delta1==1])))
	
	log2 <- sum(log(exp(-subu[delta2==1] * subx[delta2==1])
				-exp(-subv[delta2==1] * subx[delta2==1])))
				
	log3 <- -sum(subv[delta3==1] * subx[delta3==1])
	
	res <- -log1 - log2 - log3
	
	res
}
########################################################################

################ derivative of the log likelihood function  ##############################################

gradient <- function(eta) 
{
    
	eta <- matrix(c(eta), (n.total),1)
	alpha<-eta[(1):(n.spline),]
	subu<-c(t(exp(alpha))%*%ispu)
	subv<-c(t(exp(alpha))%*%ispv)
	
	if (n.total > n.spline)	{
		beta<-eta[(n.spline+1):(n.total),]
		subx<-exp(c(t(beta) %*% z + fixsubx))	
	} else
		subx <- exp(c(fixsubx))
		
	diag_alpha<-diag(c(exp(alpha)),n.spline)
	
		
	### derivative for alpha
	ggalpha11 <- t(diag_alpha %*% ispu[,delta1==1]) * 
	exp(-subu[delta1==1] * subx[delta1==1]) * subx[delta1==1]
	
	ggalpha12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
	
	ggalpha1 <- number.ones1 %*% (ggalpha11 / ggalpha12)
		
	ggalpha21 <- -t(diag_alpha %*% ispu[,delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) * subx[delta2==1]
	
	ggalpha22 <- t(diag_alpha %*% ispv[,delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) * subx[delta2==1]
	
	ggalpha23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * subx[delta2==1]) 
	
	ggalpha2 <- number.ones2 %*% ((ggalpha21 + ggalpha22) / ggalpha23)
	
	ggalpha3 <- -number.ones3 %*% (t(diag_alpha %*% ispv[,delta3==1]) * subx[delta3==1])
	
	ggalpha <- ggalpha1 + ggalpha2 + ggalpha3
	
	### derivarive for beta
	
	if (n.total > n.spline) {
		ggbeta11 <- t(z[,delta1==1]) * exp(-subu[delta1==1] * subx[delta1==1]) * 
		subu[delta1==1] * subx[delta1==1]
		
		ggbeta12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
		
		ggbeta1 <- number.ones1 %*% (ggbeta11 / ggbeta12)
		
		ggbeta21 <- -t(z[,delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) * 
		subu[delta2==1] * subx[delta2==1]
		
		ggbeta22 <- t(z[,delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) * 
		subv[delta2==1] * subx[delta2==1]
		
		ggbeta23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * subx[delta2==1])
		
		ggbeta2 <- number.ones2 %*% ((ggbeta21 + ggbeta22) / ggbeta23)
		
		ggbeta3 <- -number.ones3 %*% (t(z[,delta3==1]) * subv[delta3==1] * subx[delta3==1])
		
		ggbeta <- ggbeta1 + ggbeta2 + ggbeta3
		
		res <- -c(c(ggalpha), c(ggbeta))
	} else
		res <- -c(ggalpha)
		
	res
		
}
##############################################################################################################

##################### for variance estimation based on Zhang et al. (2010) ##############################
observ.beta <- function(eta)
{
	eta <- matrix(c(eta), (n.variable+n.spline),1)
	alpha<-eta[(1):(n.spline),]
	beta<-eta[(n.spline+1):(n.variable+n.spline),]
	
	subu<-c(t(exp(alpha))%*%ispu)
	subv<-c(t(exp(alpha))%*%ispv)
	subx<-exp(c(t(beta)%*%z))	
	
	#part 1 beta observation
	
	obeta11 <- t(z[,delta1==1]) * exp(-subu[delta1==1] * subx[delta1==1]) * 
	subu[delta1==1] * subx[delta1==1]
	
	obeta12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
	
	obeta1 <- t((obeta11 / obeta12)) %*% (obeta11 / obeta12)
	
	#part 2 beta observation
		
	obeta21 <- -t(z[,delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) * 
	subu[delta2==1] * subx[delta2==1]
	
	obeta22 <- t(z[,delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) * 
	subv[delta2==1] * subx[delta2==1]
	
	obeta23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * subx[delta2==1])
	
	obeta2 <- t((obeta21 + obeta22) / obeta23) %*% ((obeta21 + obeta22) / obeta23)	
	
	# part 3 beta observation
	
	obeta31 <- -t(z[,delta3==1]) * subv[delta3==1] * subx[delta3==1]
	
	obeta3 <- t(obeta31) %*% obeta31
	
	# total beta observation
	
	res <- (obeta1 + obeta2 + obeta3) / size
	
	res
	
}


observ.beta.lambda <- function(eta)
{
	
	eta <- matrix(c(eta), (n.variable+n.spline),1)
	alpha<-eta[(1):(n.spline),]
	beta<-eta[(n.spline+1):(n.variable+n.spline),]
	
	subu<-c(t(exp(alpha))%*%ispu)
	subv<-c(t(exp(alpha))%*%ispv)
	subx<-exp(c(t(beta)%*%z))	
	
	#part 1 beta cross Lambda observation
	
	obeta11 <- t(z[,delta1==1]) * exp(-subu[delta1==1] * subx[delta1==1]) * 
	subu[delta1==1] * subx[delta1==1]
	
	obeta12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
	
	olambda11 <- t(bspu[, delta1==1]) * exp(-subu[delta1==1] * subx[delta1==1]) *
	subx[delta1==1]
	
	olambda12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
	
	obeta.lambda1 <- t((obeta11 / obeta12)) %*% (olambda11 / olambda12)

	#part 2 beta cross Lambda observation

	obeta21 <- -t(z[,delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) * 
	subu[delta2==1] * subx[delta2==1]
	
	obeta22 <- t(z[,delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) * 
	subv[delta2==1] * subx[delta2==1]
	
	obeta23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * 
	subx[delta2==1])
	
	olambda21 <- -t(bspu[, delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) *
	subx[delta2==1]
	
	olambda22 <- t(bspv[, delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) *
	subx[delta2==1]
	
	olambda23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * 
	subx[delta2==1])
	
	obeta.lambda2 <- t((obeta21 + obeta22) / obeta23) %*% ((olambda21 + olambda22) / olambda23)
	
	# part 3 beta cross Lambda observation
	
	obeta31 <- -t(z[,delta3==1]) * subv[delta3==1] * subx[delta3==1]
	
	olambda31 <- -t(bspv[,delta3==1]) * subx[delta3==1]
	
	obeta.lambda3 <- t(obeta31) %*% olambda31
	
	# total beta cross Lambda observation
	
	res <- (obeta.lambda1 + obeta.lambda2 + obeta.lambda3) / size
	
	res
	
}		

observ.lambda <- function(eta)
{
	
	eta <- matrix(c(eta), (n.variable+n.spline),1)
	alpha<-eta[(1):(n.spline),]
	beta<-eta[(n.spline+1):(n.variable+n.spline),]
	
	subu<-c(t(exp(alpha))%*%ispu)
	subv<-c(t(exp(alpha))%*%ispv)
	subx<-exp(c(t(beta)%*%z))	
	
	#part 1  Lambda observation
	
	olambda11 <- t(bspu[, delta1==1]) * exp(-subu[delta1==1] * subx[delta1==1]) *
	subx[delta1==1]
	
	olambda12 <- 1 - exp(-subu[delta1==1] * subx[delta1==1])
	
	olambda1 <- t((olambda11 / olambda12)) %*% (olambda11 / olambda12)

	#part 2  Lambda observation

	olambda21 <- -t(bspu[, delta2==1]) * exp(-subu[delta2==1] * subx[delta2==1]) *
	subx[delta2==1]
	
	olambda22 <- t(bspv[, delta2==1]) * exp(-subv[delta2==1] * subx[delta2==1]) *
	subx[delta2==1]
	
	olambda23 <- exp(-subu[delta2==1] * subx[delta2==1]) - exp(-subv[delta2==1] * 
	subx[delta2==1])
	
	olambda2 <- t((olambda21 + olambda22) / olambda23) %*% ((olambda21 + olambda22) / olambda23)
	
	# part 3 Lambda observation
	
	olambda31 <- -t(bspv[,delta3==1]) * subx[delta3==1]
	
	olambda3 <- t(olambda31) %*% olambda31
	
	# total beta observation
	
	res <- (olambda1 + olambda2 + olambda3) / size
	
	res
	
}	
########################################################################################################	

################### knot sequence of spline functions ####################################
get_knots <- function(ctu, ctv, delta1, delta2, delta3)
{
	knot_pre<-c(ctu[delta1==1], ctu[delta1==1], ctu[delta2==1], ctv[delta2==1], ctv[delta3==1], ctv[delta3==1])

	### qurtile knots
	qmin <- 0
	q1 <- quantile(knot_pre, (1 : (n.spline + 1 - spline.ord))/(n.spline + 2 - spline.ord))
	qmax <- max(knot_pre) + 0.001
	knot <- c(rep(qmin, spline.ord), q1, rep(qmax, spline.ord))

	knot
}
############################################################################################

############# EM algorithm function reduced from a function in R library ICsurv ########################
ICsurv.EM <- function(d1, d2, d3, Li, Ri, Xp, n.int, order, g0, b0, tol, t.seq, equal = TRUE) 
{
  P <- length(b0)
  L <- length(g0)
  N <- length(d1)
  Li[d1 == 1] <- Ri[d1 == 1]
  Ri[d3 == 1] <- Li[d3 == 1]
  ti <- c(Li[d1 == 0], Ri[d3 == 0])
  if (equal == TRUE) {
    ti.max <- max(ti) + 1e-05
    ti.min <- min(ti) - 1e-05
    knots <- seq(ti.min, ti.max, length.out = (n.int + 2))
  }
  if (equal == FALSE) {
    id <- seq(0, 1, length.out = (n.int + 2))
    id <- id[-c(1, (n.int + 2))]
    ti.max <- max(ti) + 1e-05
    ti.min <- min(ti) - 1e-05
    knots <- c(ti.min, quantile(ti, id), ti.max)
  }
  bRi <- t(Ispline(x = Ri, order = order, knots = knots))
  bLi <- t(Ispline(x = Li, order = order, knots = knots))
  bt <- t(Ispline(x = t.seq, order = order, knots = knots))
  Q1 <- function(b0, b1, g0, Xp, bRi, bLi, d1, d2, d3, L) {
    g1 <- rep(-99, L)
    xb0 <- Xp %*% b0
    xb1 <- Xp %*% b1
    dz <- 1 - exp(-(bRi %*% g0) * exp(xb0))
    dw <- 1 - exp(-(bRi %*% g0 - bLi %*% g0) * exp(xb0))
    dw[d2 == 0] = 1
    EZil <- t(t(d1 * bRi) * g0) * as.vector(exp(xb0)/dz)
    EWil <- t(t(d2 * (bRi - bLi)) * g0) * as.vector(exp(xb0)/dw)
    num <- EZil + (d2 + d3) * EWil
    den <- ((d1 + d2) * bRi + d3 * bLi) * as.vector(exp(xb1))
    g1 <- apply(num, 2, sum)/apply(den, 2, sum)
    return(g1)
  }
  Q2 <- function(b1, b0, g0, Xp, bRi, bLi, d1, d2, d3, L) {
    xb0 <- Xp %*% b0
    xb1 <- Xp %*% b1
    dz <- 1 - exp(-(bRi %*% g0) * exp(xb0))
    dw <- 1 - exp(-(bRi %*% g0 - bLi %*% g0) * exp(xb0))
    dw[d2 == 0] = 1
    EZi <- d1 * (bRi %*% g0) * exp(xb0)/dz
    EWi <- d2 * (bRi %*% g0 - bLi %*% g0) * exp(Xp %*% b0)/dw
    p1 <- sum((EZi + EWi) * (Xp %*% b1))
    EZil <- t(t(d1 * bRi) * g0) * as.vector(exp(xb0)/dz)
    EWil <- t(t(d2 * (bRi - bLi)) * g0) * as.vector(exp(xb0)/dw)
    num <- EZil + (d2 + d3) * EWil
    den <- ((d1 + d2) * bRi + d3 * bLi) * as.vector(exp(xb1))
    g1 <- apply(num, 2, sum)/apply(den, 2, sum)
    p2 <- sum(t(EZil) * log(g1) + t(EWil) * log(g1))
    p3 <- sum(((d1 + d2) * (bRi %*% g1) + d3 * (bLi %*% g1)) * 
                exp(xb1))
    res <- -(p1 + p2 - p3)
    return(res)
  }
  b1 <- optim(b0, Q2, method = "Nelder-Mead", b0 = b0, g0 = g0, 
              Xp = Xp, bRi = bRi, bLi = bLi, d1 = d1, d2 = d2, d3 = d3, 
              L = L)$par
  g1 <- Q1(b0, b1, g0, Xp, bRi, bLi, d1, d2, d3, L)
  while (max(abs(c(b0, g0) - c(b1, g1))) > tol) {
    b0 <- b1
    g0 <- g1
    b1 <- optim(b0, Q2, method = "Nelder-Mead", b0 = b0, 
                g0 = g0, Xp = Xp, bRi = bRi, bLi = bLi, d1 = d1, 
                d2 = d2, d3 = d3, L = L)$par
    g1 <- Q1(b0, b1, g0, Xp, bRi, bLi, d1, d2, d3, L)
  }
  hz <- bt %*% g1
  return(list(b = b1, g = g1, hz = hz))
}
###########################################################################################

################### variance estimation functions in R library ICsurv ########################
fast.PH.Louis.ICsurv <- function(b, g, bLi, bRi, d1, d2, d3, Xp) 
{
  xb <- Xp %*% b
  A <- -t(as.vector(((d1 + d2) * bRi %*% g + d3 * bLi %*% g) * 
                      exp(xb)) * Xp) %*% Xp
  B <- -t(((d1 + d2) * bRi + d3 * bLi) * as.vector(exp(xb))) %*% 
    Xp
  C <- matrix(0, length(g), length(g))
  D <- rbind(cbind(A, t(B)), cbind(B, C))
  ci <- 1 - exp(-(bRi %*% g) * exp(xb))
  di <- 1 - exp(-(bRi %*% g - bLi %*% g) * exp(xb))
  di[d2 == 0] = 1
  hi <- d1 * (bRi %*% g) * exp(xb)
  ti <- d2 * (bRi %*% g - bLi %*% g) * exp(xb)
  VZi <- (hi^2/ci) * (1 - 1/ci) + hi/ci
  VWi <- (ti^2/di) * (1 - 1/di) + ti/di
  E <- t(as.vector((VZi + (d2 + d3) * VWi)) * Xp) %*% Xp
  zpil <- t(t(bRi) * g)/as.vector(bRi %*% g)
  num <- (bRi %*% g - bLi %*% g)
  num[d2 == 0] <- 1
  wpil <- t(t(bRi - bLi) * g)/as.vector(num)
  CovZilZi <- zpil * as.vector((hi/ci) * (1 + hi - hi/ci))
  CovWilWi <- wpil * as.vector((ti/di) * (1 + ti - ti/di))
  F <- t(CovZilZi + (d2 + d3) * CovWilWi) %*% Xp/g
  Covz <- t(zpil * as.vector((hi^2/ci) * (1 - 1/ci))) %*% zpil
  Covw <- t(wpil * as.vector((d2 + d3) * (ti^2/di) * (1 - 1/di))) %*% 
    wpil
  G <- t((Covz + Covw)/g)/g
  H <- rbind(cbind(E, t(F)), cbind(F, G))
  hess <- -(D + H)
  return(hess)
}

PH.Louis.ICsurv <- function(b, g, bLi, bRi, d1, d2, d3, Xp) 
{
  xb <- Xp %*% b
  A <- matrix(0, length(b), length(b))
  B <- matrix(0, length(g), length(b))
  C <- matrix(0, length(g), length(g))
  for (i in 1:length(b)) {
    for (j in 1:length(b)) {
      A[i, j] <- -sum(((d1 + d2) * bRi %*% g + d3 * bLi %*% 
                         g) * exp(xb) * (Xp[, i] * Xp[, j]))
    }
  }
  for (i in 1:length(b)) {
    for (j in 1:length(g)) {
      B[j, i] <- -sum(((d1 + d2) * bRi[, j] + d3 * bLi[, 
                                                       j]) * exp(xb) * Xp[, i])
    }
  }
  ci <- 1 - exp(-(bRi %*% g) * exp(xb))
  di <- 1 - exp(-(bRi %*% g - bLi %*% g) * exp(xb))
  di[d2 == 0] = 1
  hi <- d1 * (bRi %*% g) * exp(xb)
  ti <- d2 * (bRi %*% g - bLi %*% g) * exp(xb)
  D <- rbind(cbind(A, t(B)), cbind(B, C))
  E <- matrix(0, length(b), length(b))
  F <- matrix(0, length(g), length(b))
  G <- matrix(0, length(g), length(g))
  VZi <- (hi^2/ci) * (1 - 1/ci) + hi/ci
  VWi <- (ti^2/di) * (1 - 1/di) + ti/di
  for (i in 1:length(b)) {
    for (j in 1:length(b)) {
      E[i, j] <- sum((VZi + (d2 + d3) * VWi) * Xp[, i] * 
                       Xp[, j])
    }
  }
  for (i in 1:length(b)) {
    for (j in 1:length(g)) {
      zpil <- bRi[, j] * g[j]/(bRi %*% g)
      num <- (bRi %*% g - bLi %*% g)
      num[d2 == 0] <- 1
      wpil <- (bRi[, j] - bLi[, j]) * g[j]/num
      CovZilZi <- zpil * (hi/ci) * (1 + hi - hi/ci)
      CovWilWi <- wpil * (ti/di) * (1 + ti - ti/di)
      F[j, i] <- 1/(g[j]) * sum((CovZilZi + (d2 + d3) * 
                                   CovWilWi) * Xp[, i])
    }
  }
  for (i in 1:length(g)) {
    for (j in 1:length(g)) {
      zpi <- bRi[, i] * g[i]/(bRi %*% g)
      num <- (bRi %*% g - bLi %*% g)
      num[d2 == 0] <- 1
      wpi <- (bRi[, i] - bLi[, i]) * g[i]/num
      zpj <- bRi[, j] * g[j]/(bRi %*% g)
      num <- (bRi %*% g - bLi %*% g)
      num[d2 == 0] <- 1
      wpj <- (bRi[, j] - bLi[, j]) * g[j]/num
      Covz <- zpi * zpj * (hi^2/ci) * (1 - 1/ci)
      Covw <- wpi * wpj * (ti^2/di) * (1 - 1/di)
      G[i, j] <- 1/(g[i] * g[j]) * sum(Covz + (d2 + d3) * 
                                         Covw)
    }
  }
  H <- rbind(cbind(E, t(F)), cbind(F, G))
  hess <- -(D + H)
  return(hess)
}
##########################################################################################
