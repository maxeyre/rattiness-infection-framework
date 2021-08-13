

rattiness.epi.predict.alt <- function(control, 
                              rat, human, par_hat,
                              n.sim, burnin, thin,
                              grid.pred, crs.val,
                              get.raster, get.values,
                              n.pred.samples,
                              pred.int=c(0.025,0.975),
                              pred.human.N = TRUE){
  
  sp <- function(x,k) max(0,x-k)
  sp <- Vectorize(sp)
  
  inv.logit = function (x) exp(x)/(1 + exp(x))
  
  # interpret xi controls
  multi.xi.on <- control$multi.xi.on[1]
  xi.var <- control$xi.var[1]
  
  # interpret controls for random effects
  with_Ui <- control$with_Ui[1] # Rattiness nugget effect (TRUE or FALSE)
  with_human_S <- control$with_human_S[1] # Human-side S(x)
  with_human_N <- control$with_human_N[1] # Human-side nugget effect
  
  if(with_human_S ==TRUE){
    stop("Error: prediction function not yet written for Gaussian Process in human-side of framework")
  }
  
  if(get.raster == FALSE & get.values == FALSE){
    stop("Error: one of 'get.raster' or 'get.values' must be == TRUE or there will be no output")
  }
  
  # interpret controls for model covariates and random effects
  cov.rat <- control$rat[!is.na(control$rat)]
  cov.rat.sub <- cov.rat[str_detect(cov.rat, "sp.", negate=TRUE) & str_detect(cov.rat, "as.factor", negate=TRUE)]
  cov.rat.sub <- unique(c(cov.rat.sub, unlist(str_extract_all(cov.rat[str_detect(cov.rat, ",", negate=TRUE)],  "(?<=\\().+?(?=\\))"))))
  
  cov.human <- control$human[!is.na(control$human)]
  if(multi.xi.on == TRUE){cov.human <- c(cov.human, paste0("as.factor(",xi.var,")"))}
  cov.human.sub <- cov.human[str_detect(cov.human, "sp.", negate=TRUE) & str_detect(cov.human, "as.factor", negate=TRUE)]
  cov.human.sub <- unique(c(cov.human.sub, unlist(str_extract_all(cov.human[str_detect(cov.human, ",", negate=TRUE)],  "(?<=\\().+?(?=\\))")), cov.rat.sub))
  
  formula.rat <- as.formula(paste("~-1+",paste(cov.rat,collapse="+")))
  formula.human <- as.formula(paste("outcome ~",paste(cov.human,collapse="+")))
  formula.human.pred <- as.formula(paste("~",paste(cov.human,collapse="+")))
  
  # restrict rat and human datasets to only observations without NAs
  rat <- na.omit(rat[,c("X","Y","data_type","outcome","offset","offset_req", cov.rat.sub)])
  human <- na.omit(human[,c("X","Y","outcome",cov.human.sub)])
  
  # rattiness covariate model matrix at all locations
  rat.human <- bind_rows(rat[,cov.rat.sub], human[,cov.rat.sub])
  D.aux <- as.matrix(model.matrix(formula.rat, data=rat.human))
  
  ## coordinates and distances
  coords.rat <- as.matrix(rat[,c("X","Y")])
  coords.human <- as.matrix(human[,c("X","Y")])
  coords.set <- data.frame(rbind(coords.rat,coords.human))
  colnames(coords.set)  <- c("X","Y")
  id.rat <- 1:nrow(coords.rat)
  id.human <- (nrow(coords.rat)+1):(nrow(coords.rat)+nrow(coords.human))
  ID <- create.ID.coords(coords.set,~X+Y)
  coords <- unique(coords.set)
  U <- as.matrix(dist(coords))
  
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
  
  ## find all unique ID values in human locations
  loc.in.human <- which(sapply(1:N,function(x) any(ID[id.human]==x)))
  
  # human
  glm.fit <- glm(formula.human, data=human,family=binomial,x=TRUE)
  D.human.us <- glm.fit$x 
  p.human <- ncol(D.human.us)
  # standardise human covariate values
  D.h <- scale(D.human.us[,2:p.human])
  D.human <- cbind(D.human.us[,1], D.h)
  
  mean.D.h <- c(0,attr(D.h, "scaled:center")) # 0 and 1 are just to not standardise the intercept
  sd.D.h <- c(1,attr(D.h, "scaled:scale"))
  
  ## multiple xis
  if(multi.xi.on == TRUE){
    n.xi <- length(levels(as.factor(human[,xi.var])))
    xi.levels <- as.numeric(levels(as.factor(human[,xi.var])))
  }else{n.xi <- 1; xi.levels <- 1}
  
  # find locations (out of human locations 596:1262) in each xi level
  loc.in.xi <- lapply(1:length(xi.levels), function(x) unique(ID[id.human][which(human[,xi.var]==xi.levels[x])]))
  
  # ID which people are in each xi level
  if(n.xi == 2){
    ppl.in.xi1 <- ID[id.human] %in% loc.in.xi[[1]]
    ppl.in.xi2 <- ID[id.human] %in% loc.in.xi[[2]]
  }
  if(n.xi == 3){
    ppl.in.xi1 <- ID[id.human] %in% loc.in.xi[[1]]
    ppl.in.xi2 <- ID[id.human] %in% loc.in.xi[[2]]
    ppl.in.xi3 <- ID[id.human] %in% loc.in.xi[[3]]
  }
  
  par0 <- par_hat
  
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
  gamma0 <- par0[(p+12):(p+p.human+11)] # human covariates
  xi0 <- par0[(p+p.human+12):(p+p.human+12+n.xi-1)] # coefficient for Rattiness
  
  if(with_human_S==TRUE & with_human_N==FALSE){
    omega2.0 <- exp(par0[p+p.human+13+n.xi-1])
    zeta0 <- exp(par0[p+p.human+14+n.xi-1])
    ifelse(with_Ui == TRUE, psi0 <- exp(par0[p+p.human+15+n.xi-1])/(1+exp(par0[p+p.human+15+n.xi-1])), psi0 <- 1)
  }
  
  if(with_human_S==FALSE & with_human_N==TRUE){
    omega2_nugg0 <- exp(par0[p+p.human+13+n.xi-1])
    ifelse(with_Ui == TRUE, psi0 <- exp(par0[p+p.human+14+n.xi-1])/(1+exp(par0[p+p.human+14+n.xi-1])), psi0 <- 1)
  }
  
  if(with_human_S==TRUE & with_human_N==TRUE){
    omega2.0 <- exp(par0[p+p.human+13+n.xi-1])
    zeta0 <- exp(par0[p+p.human+14+n.xi-1])
    omega2_nugg0 <- exp(par0[p+p.human+15+n.xi-1])
    ifelse(with_Ui == TRUE, psi0 <- exp(par0[p+p.human+16+n.xi-1])/(1+exp(par0[p+p.human+16+n.xi-1])), psi0 <- 1)
  }
  
  if(with_human_S==FALSE & with_human_N==FALSE){
    ifelse(with_Ui == TRUE, psi0 <- exp(par0[p+p.human+13+n.xi-1])/(1+exp(par0[p+p.human+13+n.xi-1])), psi0 <- 1)
  }
  
  Sigma0 <- cov.matrix.setup(psi0, U, phi0, with_human_N, with_human_S, 
                             loc.in.xi, loc.in.human, 
                             xi0, omega2_nugg0, omega2.0, zeta0,
                             n.xi)
  
  Sigma0.inv <- solve(Sigma0) # Inverse of this matrix
  
  mu0 <- as.numeric(D%*%beta0) # covariates*coefficients section of R
  mu0.human <- as.numeric(D.human%*%gamma0) # covariates*coefficients for human section
  y <- rat$outcome
  z <- human$outcome
  offset <- rat$offset
  
  ID1 <- sort(unique(ID[id.rat][ind1]))
  ID2 <- sort(unique(ID[id.rat][ind2]))
  ID3 <- sort(unique(ID[id.rat][ind3]))
  ID4 <- sort(unique(ID[id.rat][ind4]))
  ID5 <- sort(unique(ID[id.rat][ind5]))
  ID6 <- sort(unique(ID[id.human]))
  
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
    
    # Human
    eta.human <- mu0.human
    if(n.xi==1){
      eta.human <- eta.human + xi0*R[ID[id.human]]
    }
    if(n.xi==2){
      eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*R[ID[id.human][ppl.in.xi1]] 
      eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*R[ID[id.human][ppl.in.xi2]]
    }
    if(n.xi==3){
      eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*R[ID[id.human][ppl.in.xi1]] 
      eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*R[ID[id.human][ppl.in.xi2]]
      eta.human[ppl.in.xi3] <- eta.human[ppl.in.xi3] + xi0[3]*R[ID[id.human][ppl.in.xi3]]
    }
    
    llik.human <- sum(z*eta.human-log(1+exp(eta.human)))
    
    diff.R <- R-mu0
    
    out <- as.numeric(-0.5*t(diff.R)%*%Sigma0.inv%*%diff.R)+
      llik1+llik2+llik3+llik4+llik5+llik.human
    as.numeric(out)
  }
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
    
    # Human
    eta.human <- mu0.human
    if(n.xi==1){
      eta.human <- eta.human + xi0*R[ID[id.human]]
      der.tot[ID6] <- der.tot[ID6]+
        tapply((z-exp(eta.human)/(1+exp(eta.human)))*xi0,ID[id.human],sum)
    }
    if(n.xi==2){
      eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*R[ID[id.human][ppl.in.xi1]] 
      eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*R[ID[id.human][ppl.in.xi2]]
      
      der.tot[loc.in.xi[[1]]] <- der.tot[loc.in.xi[[1]]]+
        tapply((z[ppl.in.xi1]-exp(eta.human[ppl.in.xi1])/(1+exp(eta.human[ppl.in.xi1])))*xi0[1],ID[id.human][ppl.in.xi1],sum)
      der.tot[loc.in.xi[[2]]] <- der.tot[loc.in.xi[[2]]]+
        tapply((z[ppl.in.xi2]-exp(eta.human[ppl.in.xi2])/(1+exp(eta.human[ppl.in.xi2])))*xi0[2],ID[id.human][ppl.in.xi2],sum)
    }
    if(n.xi==3){
      eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*R[ID[id.human][ppl.in.xi1]] 
      eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*R[ID[id.human][ppl.in.xi2]]
      eta.human[ppl.in.xi3] <- eta.human[ppl.in.xi3] + xi0[3]*R[ID[id.human][ppl.in.xi3]]
      
      der.tot[loc.in.xi[[1]]] <- der.tot[loc.in.xi[[1]]]+
        tapply((z[ppl.in.xi1]-exp(eta.human[ppl.in.xi1])/(1+exp(eta.human[ppl.in.xi1])))*xi0[1],ID[id.human][ppl.in.xi1],sum)
      der.tot[loc.in.xi[[2]]] <- der.tot[loc.in.xi[[2]]]+
        tapply((z[ppl.in.xi2]-exp(eta.human[ppl.in.xi2])/(1+exp(eta.human[ppl.in.xi2])))*xi0[2],ID[id.human][ppl.in.xi2],sum)
      der.tot[loc.in.xi[[3]]] <- der.tot[loc.in.xi[[3]]]+
        tapply((z[ppl.in.xi3]-exp(eta.human[ppl.in.xi3])/(1+exp(eta.human[ppl.in.xi3])))*xi0[3],ID[id.human][ppl.in.xi3],sum)
    }
    
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
    
    # Human
    eta.human <- mu0.human
    
    if(n.xi==1){
      eta.human <- eta.human + xi0*R[ID[id.human]]
      hess.tot[ID6] <- hess.tot[ID6]+
        tapply(-(exp(eta.human)/((1+exp(eta.human))^2))*xi0^2,ID[id.human],sum)
    }
    if(n.xi==2){
      eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*R[ID[id.human][ppl.in.xi1]] 
      eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*R[ID[id.human][ppl.in.xi2]]
      
      hess.tot[loc.in.xi[[1]]] <- hess.tot[loc.in.xi[[1]]]+
        tapply(-(exp(eta.human[ppl.in.xi1])/((1+exp(eta.human[ppl.in.xi1]))^2))*xi0[1]^2,ID[id.human][ppl.in.xi1],sum)
      hess.tot[loc.in.xi[[2]]] <- hess.tot[loc.in.xi[[2]]]+
        tapply(-(exp(eta.human[ppl.in.xi2])/((1+exp(eta.human[ppl.in.xi2]))^2))*xi0[2]^2,ID[id.human][ppl.in.xi2],sum)
    }
    if(n.xi==3){
      eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*R[ID[id.human][ppl.in.xi1]] 
      eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*R[ID[id.human][ppl.in.xi2]]
      eta.human[ppl.in.xi3] <- eta.human[ppl.in.xi3] + xi0[3]*R[ID[id.human][ppl.in.xi3]]
      
      hess.tot[loc.in.xi[[1]]] <- hess.tot[loc.in.xi[[1]]]+
        tapply(-(exp(eta.human[ppl.in.xi1])/((1+exp(eta.human[ppl.in.xi1]))^2))*xi0[1]^2,ID[id.human][ppl.in.xi1],sum)
      hess.tot[loc.in.xi[[2]]] <- hess.tot[loc.in.xi[[2]]]+
        tapply(-(exp(eta.human[ppl.in.xi2])/((1+exp(eta.human[ppl.in.xi2]))^2))*xi0[2]^2,ID[id.human][ppl.in.xi2],sum)
      hess.tot[loc.in.xi[[3]]] <- hess.tot[loc.in.xi[[3]]]+
        tapply(-(exp(eta.human[ppl.in.xi3])/((1+exp(eta.human[ppl.in.xi3]))^2))*xi0[3]^2,ID[id.human][ppl.in.xi3],sum)
    }
    
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
    
    # Human
    eta.human <- mu0.human
    if(n.xi==1){
      eta.human <- eta.human + xi0*R[ID[id.human]]
      der.tot[ID6] <- der.tot[ID6]+
        tapply((z-exp(eta.human)/(1+exp(eta.human)))*xi0,ID[id.human],sum)
    }
    if(n.xi==2){
      eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*R[ID[id.human][ppl.in.xi1]] 
      eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*R[ID[id.human][ppl.in.xi2]]
      
      der.tot[loc.in.xi[[1]]] <- der.tot[loc.in.xi[[1]]]+
        tapply((z[ppl.in.xi1]-exp(eta.human[ppl.in.xi1])/(1+exp(eta.human[ppl.in.xi1])))*xi0[1],ID[id.human][ppl.in.xi1],sum)
      der.tot[loc.in.xi[[2]]] <- der.tot[loc.in.xi[[2]]]+
        tapply((z[ppl.in.xi2]-exp(eta.human[ppl.in.xi2])/(1+exp(eta.human[ppl.in.xi2])))*xi0[2],ID[id.human][ppl.in.xi2],sum)
    }
    if(n.xi==3){
      eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*R[ID[id.human][ppl.in.xi1]] 
      eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*R[ID[id.human][ppl.in.xi2]]
      eta.human[ppl.in.xi3] <- eta.human[ppl.in.xi3] + xi0[3]*R[ID[id.human][ppl.in.xi3]]
      
      der.tot[loc.in.xi[[1]]] <- der.tot[loc.in.xi[[1]]]+
        tapply((z[ppl.in.xi1]-exp(eta.human[ppl.in.xi1])/(1+exp(eta.human[ppl.in.xi1])))*xi0[1],ID[id.human][ppl.in.xi1],sum)
      der.tot[loc.in.xi[[2]]] <- der.tot[loc.in.xi[[2]]]+
        tapply((z[ppl.in.xi2]-exp(eta.human[ppl.in.xi2])/(1+exp(eta.human[ppl.in.xi2])))*xi0[2],ID[id.human][ppl.in.xi2],sum)
      der.tot[loc.in.xi[[3]]] <- der.tot[loc.in.xi[[3]]]+
        tapply((z[ppl.in.xi3]-exp(eta.human[ppl.in.xi3])/(1+exp(eta.human[ppl.in.xi3])))*xi0[3],ID[id.human][ppl.in.xi3],sum)
    }
    
    diff.W <- W-mu.W
    
    as.numeric(-Sigma.W.inv%*%diff.W+
                 t(Sigma.sroot)%*%der.tot)
  }
  
  # MCMC samples
  langevin.MCMC.samples <- function(N, n.sim, n.samples, 
                                    estim, burnin, thin, Sigma.sroot){
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
      W.prop <- mean.curr+h*rnorm(N) 
      R.prop <-  as.numeric(Sigma.sroot%*%W.prop+estim$par) # transform to R scale
      mean.prop <- as.numeric(W.prop + (h^2/2)*lang.grad(W.prop,R.prop)) 
      lp.prop <- cond.dens.W(W.prop,R.prop) # density at proposal value
      
      dprop.curr <- -sum((W.prop-mean.curr)^2)/(2*(h^2)) 
      dprop.prop <- -sum((W.curr-mean.prop)^2)/(2*(h^2))
      
      log.prob <- lp.prop+dprop.prop-lp.curr-dprop.curr 
      
      if(log(runif(1)) < log.prob) { 
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
    return(sim)
  }
  
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
    
    # Humans
    eta.human <- mu0.human
    
    if(n.xi==1){
      eta.human <- eta.human + xi0*R[ID[id.human]]
    }
    if(n.xi==2){
      eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*R[ID[id.human][ppl.in.xi1]] 
      eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*R[ID[id.human][ppl.in.xi2]]
    }
    if(n.xi==3){
      eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*R[ID[id.human][ppl.in.xi1]] 
      eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*R[ID[id.human][ppl.in.xi2]]
      eta.human[ppl.in.xi3] <- eta.human[ppl.in.xi3] + xi0[3]*R[ID[id.human][ppl.in.xi3]]
    }
    
    llik.human <- sum(z*eta.human-log(1+exp(eta.human)))
    
    diff.W <- W-mu.W
    -0.5*as.numeric(t(diff.W)%*%Sigma.W.inv%*%diff.W)+
      llik1+llik2+llik3+llik4+llik5+llik.human
  }
  
  n.samples <- (n.sim-burnin)/thin
  sim <- langevin.MCMC.samples(N, n.sim, n.samples, 
                               estim, burnin, thin, Sigma.sroot)
  
  
    ## Set up for predictions ----
  D.pred <- grid.pred
  grid.pred <- grid.pred[,c("X","Y")]
  n.pred <- nrow(grid.pred)
  U.pred <- as.matrix(pdist(grid.pred,coords))

  # rat predictor
  D.pred.rat <- as.matrix(model.matrix(formula.rat, data=D.pred))
  D.pred.rat <- sapply(1:ncol(D.pred.rat),function(i) (D.pred.rat[,i]-mean.D[i])/sd.D[i])
  mu.pred.rat <- as.numeric(D.pred.rat%*%beta0)
  
  # rat S_R(x) zero predictor
  D.pred.rat.S <- matrix(0,ncol=p,nrow=n.pred)
  mu.pred.rat.S <- as.numeric(D.pred.rat.S%*%beta0)
  
  # human predictor
  D.pred.human <- as.matrix(model.matrix(formula.human.pred, data=D.pred))
  D.pred.human <- sapply(1:ncol(D.pred.human),function(i) (D.pred.human[,i]-mean.D.h[i])/sd.D.h[i])
  mu.pred.inc <- as.numeric(D.pred.human%*%gamma0)
  
  # human S_H(x) zero predictor
  D.pred.human.S <- matrix(0, ncol = p.human, nrow = n.pred)
  mu.pred.inc.S <- as.numeric(D.pred.human.S%*%gamma0)

  ## Predictions ----
  
  # Rattiness - predict R(x)
  C <- psi0*exp(-U.pred/phi0)
  A <- C%*%Sigma0.inv
  R.pred.cond.mean <- sapply(1:n.samples,function(i) mu.pred.rat+A%*%(sim[i,]-mu0))
  R.pred.hat <- apply(R.pred.cond.mean,1,mean)
  R.pred.sd <- sqrt(psi0-apply(A*C,1,sum)) 
  R.pred.cond.mean <- NULL # remove to save memory
  # correct sd for any prediction locations which are at observed location (only needed if no nugget effect; because correlation = 1)
  if(length(which(is.na(R.pred.sd)))>0){R.pred.sd[which(is.na(R.pred.sd))] <- 0}
  
  # conditional samples
  R.pred.cond.samples <- sapply(1:n.pred, function(i) R.pred.hat[i] + R.pred.sd[i]*rnorm(n.pred.samples, 0, 1))

  # Rattiness - predict S_R(x)
  S.rat.pred.cond.mean <- sapply(1:n.samples,function(i) mu.pred.rat.S+A%*%(sim[i,]-mu0))
  S.rat.pred.hat <- apply(S.rat.pred.cond.mean,1,mean)

  # find prediction locations in each xi level
  pred.loc.in.xi <- lapply(1:length(xi.levels), function(x) unique(which(D.pred[,xi.var]==xi.levels[x])))
  
  xi.vec <- rep(NA,n.pred)
  # 1 xi level
  if(n.xi==1){
    xi.vec <- xi0[1]
  }
  
  # 2 xi levels
  if(n.xi==2){
    xi.vec[pred.loc.in.xi[[1]]] <- xi0[1]
    xi.vec[pred.loc.in.xi[[2]]] <- xi0[2]
  }
  
  # 3 xi levels
  if(n.xi==3){
    xi.vec[pred.loc.in.xi[[1]]] <- xi0[1]
    xi.vec[pred.loc.in.xi[[2]]] <- xi0[2]
    xi.vec[pred.loc.in.xi[[3]]] <- xi0[3]
  }
  
  # Using samples from R|W
  ifelse(with_human_N==TRUE & pred.human.N ==TRUE,
         human.pred.cond.samples <- sapply(1:n.pred, 
                                           function(i) 1/(1+exp(-(mu.pred.inc[i] + xi.vec[i]*R.pred.cond.samples[,i] + 
                                                                    rnorm(n.pred.samples,0,sqrt(omega2_nugg0)))))),
         human.pred.cond.samples <- sapply(1:n.pred, 
                                           function(i) 1/(1+exp(-(mu.pred.inc[i] + xi.vec[i]*R.pred.cond.samples[,i])))))
  
  human.pred.hat <- apply(human.pred.cond.samples,2,mean)
  human.pred.sd <- apply(human.pred.cond.samples,2,sd)
  human.pred.hat.int <- apply(human.pred.cond.samples,2,function(x) quantile(x,c(pred.int[1],pred.int[2])))
  
  if(get.values==TRUE){
    if(with_human_S==FALSE){
      predictions.val <- list(R.pred.cond.samples = R.pred.cond.samples, R.pred.hat = R.pred.hat, R.pred.sd = R.pred.sd,
                              S.rat.pred.hat = S.rat.pred.hat,
                              human.pred.cond.samples = human.pred.cond.samples, human.pred.hat = human.pred.hat, 
                              human.pred.sd = human.pred.sd, human.pred.hat.int = human.pred.hat.int)
      
    }else if(with_human_S==TRUE){
      predictions.val <- list(R.pred.cond.samples = R.pred.cond.samples, R.pred.hat = R.pred.hat, R.pred.sd = R.pred.sd,
                              S.rat.pred.hat = S.rat.pred.hat,
                              human.pred.cond.samples = human.pred.cond.samples, human.pred.hat = human.pred.hat, 
                              human.pred.sd = human.pred.sd, human.pred.hat.int = human.pred.hat.int, 
                              S.human.pred.cond.mean = S.human.pred.cond.mean, S.human.pred.hat = S.human.pred.hat)
    }
  }
  
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

    # Human - incidence (either with or without S_H(x) depending on previous if statement)
    r.human.hat <- rasterFromXYZ(cbind(grid.pred,human.pred.hat))
    r.human.sd <- rasterFromXYZ(cbind(grid.pred,human.pred.sd))
    crs(r.human.hat) <- crs.val
    crs(r.human.sd) <- crs.val

    # Human - prediction interval
    r.human.hat.li <- rasterFromXYZ(cbind(grid.pred,human.pred.hat.int[1,]))
    r.human.hat.ui <- rasterFromXYZ(cbind(grid.pred,human.pred.hat.int[2,]))
    crs(r.human.hat.li) <- crs.val
    crs(r.human.hat.ui) <- crs.val

    if(with_human_S==TRUE){
      # Human - only S_H(x)
      r.S.human.hat <- rasterFromXYZ(cbind(grid.pred,S.human.pred.hat))
      r.S.human.sd <- rasterFromXYZ(cbind(grid.pred,S.human.pred.sd))
      crs(r.S.human.hat) <- crs.val
      crs(r.S.human.sd) <- crs.val
    }
    
    if(with_human_S==FALSE){
      predictions.raster <- list(r.R.hat = r.R.hat, r.R.sd = r.R.sd,
                                 r.S.rat.hat = r.S.rat.hat,
                                 r.human.hat = r.human.hat, r.human.sd = r.human.sd, 
                                 r.human.hat.li = r.human.hat.li, r.human.hat.ui = r.human.hat.ui)
    }else if(with_human_S==TRUE){
      predictions.raster <- list(r.R.hat = r.R.hat, r.R.sd = r.R.sd,
                                 r.S.rat.hat = r.S.rat.hat,
                                 r.human.hat = r.human.hat, r.human.sd = r.human.sd,
                                 r.human.hat.li = r.human.hat.li, r.human.hat.ui = r.human.hat.ui,
                                 r.S.human.hat = r.S.human.hat, r.S.human.sd = r.S.human.sd, )
    }
  }
  
  if(get.raster == TRUE & get.values == TRUE){
    out <- list(predictions.val, predictions.raster)
    names(out) <- c("predictions.val", "predictions.raster")
  }
  if(get.raster == TRUE & get.values == FALSE){
    out <- list(predictions.raster)
    names(out) <- c("predictions.raster")
  }
  if(get.raster == FALSE & get.values == TRUE){
    out <- list(predictions.val)
    names(out) <- c("predictions.val")
  }

  return(out)
}
