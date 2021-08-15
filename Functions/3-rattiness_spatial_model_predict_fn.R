
# 1. make predictions

rattiness.ECO.predict <- function(control, 
                              rat, par_hat,
                              n.sim, burnin, thin,
                              grid.pred, crs.val,
                              get.raster, n.pred.samples){
  
  sp <- function(x,k) max(0,x-k)
  sp <- Vectorize(sp)
  
  # interpret controls for random effects
  with_Ui <- control$with_Ui[1] # Rattiness nugget effect (TRUE or FALSE)
  
  # interpret controls for model covariates and random effects
  cov.rat <- control$rat[!is.na(control$rat)]
  cov.rat.sub <- cov.rat[str_detect(cov.rat, "sp.", negate=TRUE) & str_detect(cov.rat, "as.factor", negate=TRUE)]
  cov.rat.sub <- unique(c(cov.rat.sub, unlist(str_extract_all(cov.rat[str_detect(cov.rat, ",", negate=TRUE)],  "(?<=\\().+?(?=\\))"))))
  
  formula.rat <- as.formula(paste("~-1+",paste(cov.rat,collapse="+")))
  
  # restrict rat dataset to only observations without NAs
  rat <- na.omit(rat[,c("X","Y","data_type","outcome","offset","offset_req", cov.rat.sub)])
  
  # rattiness covariate model matrix
  D.aux <- as.matrix(model.matrix(formula.rat, data=rat))
  
  ## coordinates and distances
  ID <- create.ID.coords(rat,~X + Y)
  coords <- rat[,c("X","Y")]
  id.rat <- 1:nrow(coords)
  coords <- unique(coords)
  U <- dist(coords)
  
  
  ## standardise rattiness covariates
  p <- ncol(D.aux)
  N <- nrow(coords)
  D <- matrix(NA,nrow=N,ncol=p)
  for(i in 1:p) {
    D[,i] <- tapply(D.aux[,i],ID,max)
  }
  D.unscale <- D
  D <- scale(D)
  mean.D <- attr(D, "scaled:center")
  sd.D <- attr(D, "scaled:scale")
  
  ind1 <- which(rat$data_type=="traps")
  ind2 <- which(rat$data_type=="plates")
  ind3 <- which(rat$data_type=="burrows")
  ind4 <- which(rat$data_type=="faeces")
  ind5 <- which(rat$data_type=="trails")
  
  ## parameters estimates
  
  # rat
  par0 <- par_hat$`final estimates`$par
    
  alpha1.0 <- par0[1]
  alpha2.0 <- par0[2]
  alpha3.0 <- par0[3]
  alpha4.0 <- par0[4]
  alpha5.0 <- par0[5]
  beta0 <- par0[11:(p+10)] # rattiness covariates
  sigma1.0 <- exp(par0[6])
  sigma2.0 <- exp(par0[7])
  sigma3.0 <- exp(par0[8])
  sigma4.0 <- exp(par0[9])
  sigma5.0 <- exp(par0[10])
  phi0 <- exp(par0[p+11]) # scale of spatial correlation in spatial Gaussian process in rattiness
  
  ifelse(with_Ui == TRUE, psi0 <- exp(par0[p+12])/(1+exp(par0[p+12])), psi0 <- 1)
  
  Sigma0 <- as.matrix(psi0*exp(-U/phi0))
  diag(Sigma0) <- 1
  
  Sigma0.inv <- solve(Sigma0) # Inverse of this matrix
  
  mu0 <- as.numeric(D%*%beta0) # covariates*coefficients section of R
  
  y <- rat$outcome
  offset <- rat$offset
  
  integrand <- function(R) {
    
    # Traps
    lambda1 <- exp(alpha1.0+sigma1.0*R[ID[id.rat][ind1]])
    prob1 <- 1-exp(-offset[ind1]*lambda1)
    llik1 <- sum(y[ind1]*log(prob1/(1-prob1))+log(1-prob1))
    
    # Plates
    eta2 <- alpha2.0+sigma2.0*R[ID[id.rat][ind2]]
    llik2 <- sum(y[ind2]*eta2-offset[ind2]*log(1+exp(eta2)))
    
    # Burrows
    lambda3 <- exp(alpha3.0+sigma3.0*R[ID[id.rat][ind3]])
    llik3 <- sum(-lambda3 + y[ind3]*log(lambda3))
    
    # Faeces
    eta4 <- alpha4.0+sigma4.0*R[ID[id.rat][ind4]]
    llik4 <- sum(y[ind4]*eta4-offset[ind4]*log(1+exp(eta4)))
    
    # Trails
    eta5 <- alpha5.0+sigma5.0*R[ID[id.rat][ind5]]
    llik5 <- sum(y[ind5]*eta5-offset[ind5]*log(1+exp(eta5)))
    
    diff.R <- R-mu0 
    out <- as.numeric(-0.5*t(diff.R)%*%Sigma0.inv%*%diff.R)+
      llik1+llik2+llik3+llik4+llik5
    as.numeric(out)
  }
  
  ID1 <- sort(unique(ID[id.rat][ind1]))
  ID2 <- sort(unique(ID[id.rat][ind2]))
  ID3 <- sort(unique(ID[id.rat][ind3]))
  ID4 <- sort(unique(ID[id.rat][ind4]))
  ID5 <- sort(unique(ID[id.rat][ind5]))
  
  grad.integrand <- function(R) {
    der.tot <- rep(0,N)
    
    # Traps
    lambda1 <- exp(alpha1.0+sigma1.0*R[ID[ind1]])
    prob1 <- 1-exp(-offset[ind1]*lambda1)
    der.prob1 <- offset[ind1]*exp(-offset[ind1]*lambda1)*lambda1*sigma1.0
    der.tot[ID1] <- der.tot[ID1]+
      tapply((y[ind1]/(prob1*(1-prob1))-1/(1-prob1))*der.prob1,ID[ind1],sum)
    
    # Plates
    eta2 <- alpha2.0+sigma2.0*R[ID[ind2]]
    der.tot[ID2] <- der.tot[ID2]+
      tapply((y[ind2]-offset[ind2]*exp(eta2)/(1+exp(eta2)))*sigma2.0,ID[ind2],sum)
    
    # Burrows
    lambda3 <- exp(alpha3.0+sigma3.0*R[ID[ind3]])
    der.tot[ID3] <- der.tot[ID3] + 
      tapply(-sigma3.0*lambda3 + y[ind3]*sigma3.0, ID[ind3],sum)
    
    # Faeces
    eta4 <- alpha4.0+sigma4.0*R[ID[ind4]]
    der.tot[ID4] <- 
      tapply((y[ind4]-offset[ind4]*exp(eta4)/(1+exp(eta4)))*sigma4.0,ID[ind4],sum)
    
    # Trails
    eta5 <- alpha5.0+sigma5.0*R[ID[ind5]]
    der.tot[ID5] <- 
      tapply((y[ind5]-offset[ind5]*exp(eta5)/(1+exp(eta5)))*sigma5.0,ID[ind5],sum)
    
    diff.R <- R-mu0
    out <- -Sigma0.inv%*%diff.R+der.tot
    as.numeric(out)
  }
  
  hessian.integrand <- function(R) {
    hess.tot <- rep(0,N)
    
    # Traps
    lambda1 <- exp(alpha1.0+sigma1.0*R[ID[ind1]])
    prob1 <- 1-exp(-offset[ind1]*lambda1)
    der.prob1 <- offset[ind1]*exp(-offset[ind1]*lambda1)*lambda1*sigma1.0
    der2.prob1 <- -((offset[ind1])^2)*exp(-offset[ind1]*lambda1)*(lambda1*sigma1.0)^2+
      offset[ind1]*exp(-offset[ind1]*lambda1)*lambda1*(sigma1.0)^2
    hess.tot[ID1] <- hess.tot[ID1]+
      tapply((y[ind1]/(prob1*(1-prob1))-1/(1-prob1))*der2.prob1+
               (y[ind1]*((2*prob1-1)/((prob1*(1-prob1))^2))-1/(1-prob1)^2)*(der.prob1^2),ID[ind1],sum)
    
    # Plates
    eta2 <- alpha2.0+sigma2.0*R[ID[ind2]]
    hess.tot[ID2] <- hess.tot[ID2]+
      tapply(-(offset[ind2]*exp(eta2)/((1+exp(eta2))^2))*sigma2.0^2,ID[ind2],sum)
    
    # Burrows
    lambda3 <- exp(alpha3.0+sigma3.0*R[ID[ind3]])
    hess.tot[ID3] <- hess.tot[ID3] + 
      tapply(-sigma3.0^2*lambda3, ID[ind3],sum)
    
    # Faeces
    eta4 <- alpha4.0+sigma1.0*R[ID[ind4]]
    hess.tot[ID4] <- hess.tot[ID4]+
      tapply(-(offset[ind4]*exp(eta4)/((1+exp(eta4))^2))*sigma4.0^2,ID[ind4],sum)
    
    # Trails
    eta5 <- alpha5.0+sigma5.0*R[ID[ind5]]
    hess.tot[ID5] <- hess.tot[ID5]+
      tapply(-(offset[ind5]*exp(eta5)/((1+exp(eta5))^2))*sigma5.0^2,ID[ind5],sum)
    
    
    out <- -Sigma0.inv
    diag(out) <- diag(out)+hess.tot
    out
  }
  
  estim <- nlminb(start=rep(0,N),
                  function(x) -integrand(x),
                  function(x) -grad.integrand(x),
                  function(x) -hessian.integrand(x))
  H <- hessian.integrand(estim$par)
  
  Sigma.sroot <- t(chol(solve(-H)))
  A <- solve(Sigma.sroot)
  Sigma.W.inv <- solve(A%*%Sigma0%*%t(A))
  mu.W <- as.numeric(A%*%(mu0-estim$par))
  
  
  cond.dens.W <- function(W,R) {
    
    # Traps
    lambda1 <- exp(alpha1.0+sigma1.0*R[ID[ind1]])
    prob1 <- 1-exp(-offset[ind1]*lambda1)
    llik1 <- sum(y[ind1]*log(prob1/(1-prob1))+log(1-prob1))
    
    # Plates
    eta2 <- alpha2.0+sigma2.0*R[ID[ind2]]
    llik2 <- sum(y[ind2]*eta2-offset[ind2]*log(1+exp(eta2)))
    
    # Burrows
    lambda3 <- exp(alpha3.0+sigma3.0*R[ID[ind3]])
    llik3 <- sum(-lambda3 + y[ind3]*log(lambda3))
    
    # Faeces
    eta4 <- alpha4.0+sigma4.0*R[ID[ind4]]
    llik4 <- sum(y[ind4]*eta4-offset[ind4]*log(1+exp(eta4)))
    
    # Trails
    eta5 <- alpha5.0+sigma5.0*R[ID[ind5]]
    llik5 <- sum(y[ind5]*eta5-offset[ind5]*log(1+exp(eta5)))
    
    diff.W <- W-mu.W
    -0.5*as.numeric(t(diff.W)%*%Sigma.W.inv%*%diff.W)+
      llik1+llik2+llik3+llik4+llik5
  }
  
  # Gradient for langevin wrt W - this is part of proposal distribution for langevin MCMC
  lang.grad <- function(W,R) {
    der.tot <- rep(0,N)
    
    # Traps
    lambda1 <- exp(alpha1.0+sigma1.0*R[ID[ind1]])
    prob1 <- 1-exp(-offset[ind1]*lambda1)
    der.prob1 <- offset[ind1]*exp(-offset[ind1]*lambda1)*lambda1*sigma1.0
    der.tot[ID1] <- der.tot[ID1]+
      tapply((y[ind1]/(prob1*(1-prob1))-1/(1-prob1))*der.prob1,ID[ind1],sum)
    
    # Plates
    eta2 <- alpha2.0+sigma2.0*R[ID[ind2]]
    der.tot[ID2] <- 
      tapply((y[ind2]-offset[ind2]*exp(eta2)/(1+exp(eta2)))*sigma2.0,ID[ind2],sum)
    
    # Burrows
    lambda3 <- exp(alpha3.0+sigma3.0*R[ID[ind3]])
    der.tot[ID3] <- der.tot[ID3] + 
      tapply(-sigma3.0*lambda3 + y[ind3]*sigma3.0, ID[ind3],sum)
    
    # Faeces
    eta4 <- alpha4.0+sigma4.0*R[ID[ind4]]
    der.tot[ID4] <- der.tot[ID4]+
      tapply((y[ind4]-offset[ind4]*exp(eta4)/(1+exp(eta4)))*sigma4.0,ID[ind4],sum)
    
    # Trails
    eta5 <- alpha5.0+sigma5.0*R[ID[ind5]]
    der.tot[ID5] <- der.tot[ID5]+
      tapply((y[ind5]-offset[ind5]*exp(eta5)/(1+exp(eta5)))*sigma5.0,ID[ind5],sum)
    
    diff.W <- W-mu.W
    
    as.numeric(-Sigma.W.inv%*%diff.W+
                 t(Sigma.sroot)%*%der.tot)
  }
  
  
  h <- 1.65/(N^(1/6))
  
  c1.h <- 0.001
  c2.h <- 0.0001
  W.curr <- rep(0,N)
  
  R.curr <- as.numeric(Sigma.sroot%*%W.curr+estim$par)
  mean.curr <- as.numeric(W.curr + (h^2/2)*lang.grad(W.curr,R.curr)) ## definition of langevin distribution (pushes you up gradient to regions of higher probably density)
  lp.curr <- cond.dens.W(W.curr,R.curr) # density of the conditional distribution of W on the log scale (our target distribution)
  acc <- 0
  n.samples <- (n.sim-burnin)/thin
  sim <- matrix(NA,nrow=n.samples,ncol=N)
  
  h.vec <- rep(NA,n.sim)
  for(i in 1:n.sim) {
    W.prop <- mean.curr+h*rnorm(N) # random walk (not quite as W.prop changes)
    R.prop <-  as.numeric(Sigma.sroot%*%W.prop+estim$par) # transform to R scale
    mean.prop <- as.numeric(W.prop + (h^2/2)*lang.grad(W.prop,R.prop)) 
    lp.prop <- cond.dens.W(W.prop,R.prop) # density at proposal value
    
    dprop.curr <- -sum((W.prop-mean.curr)^2)/(2*(h^2)) # standard definition of ratio used in metropolis-hasting algorithm
    dprop.prop <- -sum((W.curr-mean.prop)^2)/(2*(h^2))
    
    log.prob <- lp.prop+dprop.prop-lp.curr-dprop.curr # difference in log probabilities (= ratio of probabilities) # look up this ratio, look up basics of metropolis-hastings
    
    if(log(runif(1)) < log.prob) { # discuss condition for accepting proposed value of W. Point here is to 'bias' towards higher prob areas?
      acc <- acc+1
      W.curr <- W.prop
      R.curr <- R.prop
      lp.curr <- lp.prop
      mean.curr <- mean.prop
    }
    
    if( i > burnin & (i-burnin)%%thin==0) {
      sim[(i-burnin)/thin,] <- R.curr
    }
    
    h.vec[i] <- h <- max(0,h + c1.h*i^(-c2.h)*(acc/i-0.57))
    #cat("Iteration",i,"out of",n.sim,"\r")
    flush.console()
  }
  
  ## Set up for predictions ----
  D.pred <- grid.pred
  grid.pred <- grid.pred[,c("X","Y")]
  n.pred <- nrow(grid.pred)
  ID.pred <- create.ID.coords(grid.pred, coords=~X+Y)
  n.pred.unq <- nrow(unique(grid.pred))
  U.pred <- as.matrix(pdist(grid.pred,coords))
  
  # rat predictor
  D.pred.rat <- as.matrix(model.matrix(formula.rat, data=D.pred))
  D.pred.rat <- sapply(1:ncol(D.pred.rat),function(i) (D.pred.rat[,i]-mean.D[i])/sd.D[i])
  mu.pred.rat <- as.numeric(D.pred.rat%*%beta0)
  
  # rat S_R(x) zero predictor
  D.pred.rat.S <- matrix(0,ncol=p,nrow=n.pred)
  mu.pred.rat.S <- as.numeric(D.pred.rat.S%*%beta0)

  ## Predictions ----
  
  # Rattiness - predict R(x)
  C <- psi0*exp(-U.pred/phi0)
  A <- C%*%Sigma0.inv
  R.pred.cond.mean <- sapply(1:n.samples,function(i) mu.pred.rat+A%*%(sim[i,]-mu0))
  R.pred.hat <- apply(R.pred.cond.mean,1,mean)
  R.pred.sd <- sqrt(psi0-apply(A*C,1,sum)) # psi0 is the variance
  # correct sd for any prediction locations which are at observed location (only needed if no nugget effect; because correlation = 1)
  if(length(which(is.na(R.pred.sd)))>0){R.pred.sd[which(is.na(R.pred.sd))] <- 0}
  # conditional samples
  R.pred.hat.unq <- unique(R.pred.hat)
  R.pred.sd.unq <- unique(R.pred.sd)
  R.pred.cond.samples <- sapply(1:n.pred.unq, function(i) R.pred.hat.unq[i] + R.pred.sd.unq[i]*rnorm(n.pred.samples, 0, 1))
  R.pred.cond.samples <- R.pred.cond.samples[,ID.pred]
  
  # Rattiness - predict S_R(x)
  S.rat.pred.cond.mean <- sapply(1:n.samples,function(i) mu.pred.rat.S+A%*%(sim[i,]-mu0))
  S.rat.pred.hat <- apply(S.rat.pred.cond.mean,1,mean)
  
  predictions.val <- list(R.pred.cond.mean = R.pred.cond.mean, R.pred.hat = R.pred.hat, R.pred.sd = R.pred.sd,
                          R.pred.cond.samples = R.pred.cond.samples,
                          S.rat.pred.cond.mean = S.rat.pred.cond.mean, S.rat.pred.hat = S.rat.pred.hat)
  
  out <- list(predictions.val)
  names(out) <- c("predictions.val")
  
  if(get.raster == TRUE){
    ## Make rasters ----
    # Rattiness - predict R(x)
    r.R.hat <- rasterFromXYZ(cbind(grid.pred,R.pred.hat))
    r.R.sd <- rasterFromXYZ(cbind(grid.pred,R.pred.sd))
    crs(r.R.hat) <- crs.val
    crs(r.R.sd) <- crs.val
    
    # Rattiness - predict S_R(x)
    r.S.rat.hat <- rasterFromXYZ(cbind(grid.pred,S.rat.pred.hat))
    crs(r.S.rat.hat) <- crs.val
    
    
    predictions.raster <- list(r.R.hat = r.R.hat, r.R.sd = r.R.sd,
                               r.S.rat.hat = r.S.rat.hat)
    out <- list(predictions.val, predictions.raster)
    names(out) <- c("predictions.val", "predictions.raster")
  }

  return(out)
}
