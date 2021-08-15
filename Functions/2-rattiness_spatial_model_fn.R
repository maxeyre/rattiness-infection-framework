rattiness.eco.model.alt <- function(control, rat,
                                euclid.norm.method, tol, n.iter.fit, rel.tol,
                                n.sim, burnin, thin, MCMC.swap, n.MCMC.swap){
  
  if(MCMC.swap ==TRUE & (length(burnin)!=2 | length(n.sim)!=2)){
    stop("Please provide two n.sim and two burnin values (in burnin vector)")
    }
  
  if(MCMC.swap ==TRUE & !is.numeric(n.MCMC.swap)){
    stop("Please designate the number of iterations after which to implement MCMC chain length change (n.MCMC.swap)")
    }
  
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

  ind1 <- which(rat$data_type=="traps")
  ind2 <- which(rat$data_type=="plates")
  ind3 <- which(rat$data_type=="burrows")
  ind4 <- which(rat$data_type=="faeces")
  ind5 <- which(rat$data_type=="trails")
  
  ## initial parameters

  # Random effect parameters
  psi0 <- control$psi0[1] # relative contribution of S(x) and Ui in rattiness
  
  # rat
  par0.rat <- control$rat_par0[!is.na(control$rat_par0)]
  
  if(is.logical(par0.rat[1]) & par0.rat[1]==FALSE){
    ifelse(with_Ui == FALSE, par0.rat <- rep(0,10+ncol(D.aux)+1), par0.rat <- c(rep(0,10+ncol(D.aux)+1)))
  }

  # full initial parameter set
  if(with_Ui == TRUE){
    par0 <- c(par0.rat,log(psi0/(1-psi0)))
  }else {
    par0 <- c(par0.rat)
  }
  
  time1 <- Sys.time()
  
  k <- 1
  
  par.matrix <- matrix(NA,nrow = length(par0), ncol=1000)
  par.matrix[,1] <- par0
  par_rel <- matrix(NA,nrow = length(par0), ncol=1000)
  euclid.norm <- c(tol + 0.1)
  ifelse(euclid.norm.method == TRUE, iter.check <- euclid.norm[k], {iter.check <- n.iter.fit; tol <- k-1})
  
  burnin.ctrl <- burnin
  n.sim.ctrl <- n.sim
  
  burnin <- burnin.ctrl[1]
  n.sim <- n.sim.ctrl[1]
  
  while(iter.check > tol) {
    
    if(MCMC.swap == TRUE & k >= n.MCMC.swap){burnin <- burnin.ctrl[2]; n.sim <- n.sim.ctrl[2]}
    
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
    
    # acf.plot <- acf(sim[,1],plot=FALSE)
    # plot(acf.plot$lag,acf.plot$acf,type="l",xlab="lag",ylab="autocorrelation",
    #      ylim=c(-0.1,1),main="Autocorrelogram of the simulated samples")
    # for(i in 2:ncol(sim)) {
    #   acf.plot <- acf(sim[,i],plot=FALSE)
    #   lines(acf.plot$lag,acf.plot$acf)
    # }
    # abline(h=0,lty="dashed",col=2)
    
    
    log.integrand <- function(R,val) {
      
      # Traps
      lambda1 <- exp(val$alpha1+val$sigma1*R[ID[ind1]])
      prob1 <- 1-exp(-offset[ind1]*lambda1)
      llik1 <- sum(y[ind1]*log(prob1/(1-prob1))+log(1-prob1))
      
      # Plates
      eta2 <- val$alpha2+val$sigma2*R[ID[ind2]]
      llik2 <- sum(y[ind2]*eta2-offset[ind2]*log(1+exp(eta2)))
      
      # Burrows
      lambda3 <- exp(val$alpha3+val$sigma3*R[ID[ind3]])
      llik3 <- sum(-lambda3 + y[ind3]*log(lambda3))
      
      # Faeces
      eta4 <- val$alpha4+val$sigma4*R[ID[ind4]]
      llik4 <- sum(y[ind4]*eta4-offset[ind4]*log(1+exp(eta4)))
      
      # Trails
      eta5 <- val$alpha5+val$sigma5*R[ID[ind5]]
      llik5 <- sum(y[ind5]*eta5-offset[ind5]*log(1+exp(eta5)))
      
      diff.R <- R-val$mu
      out <- as.numeric(-0.5*(val$log.det.Sigma+t(diff.R)%*%val$Sigma.inv%*%diff.R))+ 
        llik1+llik2+llik3+llik4+llik5
    }
    
    compute.log.f <- function(par,ldetR=NA,R.inv=NA) {
      val <- list()
      val$alpha1 <- par[1]
      val$alpha2 <- par[2]
      val$alpha3 <- par[3]
      val$alpha4 <- par[4]
      val$alpha5 <- par[5]
      
      val$sigma1 <- exp(par[6])
      val$sigma2 <- exp(par[7])
      val$sigma3 <- exp(par[8])
      val$sigma4 <- exp(par[9])
      val$sigma5 <- exp(par[10])
      
      beta <- par[11:(p+10)]
      val$mu <- as.numeric(D%*%beta)
      phi <- exp(par[p+11])
      
      ifelse(with_Ui == TRUE, psi <- exp(par[p+12])/(1+exp(par[p+12])), psi <- 1)
      
      Sigma <- as.matrix(psi*exp(-U/phi))
      diag(Sigma) <- 1
      
      val$Sigma.inv <- solve(Sigma)
      val$log.det.Sigma <- determinant(Sigma)$modulus
      
      sapply(1:(dim(sim)[1]),function(i) log.integrand(sim[i,],val))
    }
    
    
    if(with_Ui == TRUE){
      par0 <- c(alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                log(sigma1.0),log(sigma2.0),log(sigma3.0),log(sigma4.0),log(sigma5.0),
                beta0,log(phi0),log(psi0/(1-psi0)))
      names(par0) <- c("alpha_traps","alpha_plates","alpha_burrows","alpha_faeces","alpha_trails",
                       "log(sigma_traps)","log(sigma_plates)","log(sigma_burrows)","log(sigma_faeces)","log(sigma_trails)",
                       paste0("rat__",c(colnames(D.aux))),"log(phi)","log(psi/(1-psi))")
    }else{
      par0 <- c(alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                log(sigma1.0),log(sigma2.0),log(sigma3.0),log(sigma4.0),log(sigma5.0),
                beta0,log(phi0))
      names(par0) <- c("alpha_traps","alpha_plates","alpha_burrows","alpha_faeces","alpha_trails",
                       "log(sigma_traps)","log(sigma_plates)","log(sigma_burrows)","log(sigma_faeces)","log(sigma_trails)",
                       paste0("rat__",c(colnames(D.aux))),"log(phi)")
    }
    
    
    log.f.tilde <- compute.log.f(par0)
    
    MC.log.lik <- function(par) {
      log(mean(exp(compute.log.f(par)-log.f.tilde)))
    }
    
    MC.log.lik(par0)
    
    estim.par <- nlminb(par0,
                        function(x) -MC.log.lik(x),
                        control=list(trace=1, rel.tol = rel.tol))
    
    par0 <- estim.par$par
    par.matrix[,k+1] <- par0
    par_rel[,k] <- (par.matrix[,k+1] - par.matrix[,k])/par.matrix[,k]
    euclid.norm[k+1] <- sqrt(sum(par_rel[,k]^2))
    
    k <- k + 1
    
    ifelse(euclid.norm.method == TRUE, iter.check <- euclid.norm[k], {iter.check <- n.iter.fit; tol <- k-1})
    
    time2 <- Sys.time()
    
    time.diff <- difftime(time2,time1)
    
    # output parameters with names
    if(euclid.norm.method == TRUE){
      message("Iteration #", k-1, " completed. Euclid norm = ", euclid.norm[k],". Total time elapsed: ", round(time.diff,3), units(time.diff), ". MCMC simulations: ", n.sim)
      message("Parameter estimates: ")
      message(paste(names(par0), round(par0, 5), " | "))
    }else{
      message("Iteration #", k-1, " completed."," Total time elapsed: ", round(time.diff,3), units(time.diff), ". MCMC simulations: ", n.sim)
      message("Parameter estimates: ")
      message(paste(names(par0), round(par0, 5), " | "))
    }
    
  }
  
  par.matrix <- par.matrix[, colSums(is.na(par.matrix)) != nrow(par.matrix)]
  par_rel <- par_rel[, colSums(is.na(par_rel)) != nrow(par_rel)]
  
  # unstandardised regression coeffcients
  par_us_rat <- unstandardise.coeff.fn(par0[paste0("rat__",c(colnames(D.aux)))],intercept = FALSE, x = D.unscale)
  par_us <- par0
  par_us[c(paste0("rat__",c(colnames(D.aux))))] <- c(par_us_rat)
  par.list <- list(estim.par, par_us, par.matrix, par_rel, euclid.norm, time2-time1, log.f.tilde)
  names(par.list) <- c("final estimates", "rescaled regression coefficients", "iteration parameters", "relative parameter difference", "euclid norm", "total time elapsed","log-likelihood (samples)")
  
  return(par.list)  
}
