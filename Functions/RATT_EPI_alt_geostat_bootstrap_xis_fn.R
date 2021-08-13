rattiness.epi.bootstrap.alt <- function(control, rat, human, par_hat,
                                    euclid.norm.method, tol, 
                                    n.iter.fit, rel.tol, iter.max,
                                    n.sim, burnin, thin, MCMC.swap, n.MCMC.swap){
  
  if(MCMC.swap ==TRUE & (length(burnin)!=2 | length(n.sim)!=2)){
    stop("Please provide two n.sim and two burnin values (in burnin vector)")
  }
  
  if(MCMC.swap ==TRUE & !is.numeric(n.MCMC.swap)){
    stop("Please designate the number of iterations after which to implement MCMC chain length change (n.MCMC.swap)")
  }
  
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
  
  ind1 <- which(rat$data_type=="traps")
  ind2 <- which(rat$data_type=="plates")
  ind3 <- which(rat$data_type=="burrows")
  ind4 <- which(rat$data_type=="faeces")
  ind5 <- which(rat$data_type=="trails")
  
  ## find all unique ID values in human locations
  loc.in.human <- which(sapply(1:N,function(x) any(ID[id.human]==x)))
  loc.in.rat <- which(sapply(1:N,function(x) any(ID[id.rat]==x)))
  
  # human
  glm.fit <- glm(formula.human, data=human,family=binomial,x=TRUE)
  D.human.us <- glm.fit$x 
  p.human <- ncol(D.human.us)
  
  # standardise human covariate values
  D.human <- cbind(D.human.us[,1], scale(D.human.us[,2:p.human])) 
  glm.par0.s <- standardise.coeff.fn(beta = coef(glm.fit), intercept = TRUE, x = D.human.us)
  
  ## initial parameters
  # xi
  if(multi.xi.on == TRUE){
    n.xi <- length(levels(as.factor(human[,xi.var])))
    xi.levels <- as.numeric(levels(as.factor(human[,xi.var])))
  }else{n.xi <- 1; xi.levels <- 1}
  xi0 <- control$xi0[1:n.xi] # contribution of rattiness to human linear predictor
  
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
  
  # fitted parameters & make covariance matrix
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
  
  # Sampling W_R for rats and W_H for humans
  randeffects <- data.frame(ID=unique(ID), value= as.vector(rmvnorm(n=1, mean=rep(0,N), sigma=Sigma0)))
  
  mu0 <- as.numeric(D%*%beta0)
  Rj <- mu0 + randeffects$value
  
  # Create new simulated rat data
  rat_new <- rat
  rat_new$Rj <- Rj[ID[id.rat]]
  
  rat_new$alpha <-c(0)
  rat_new$sigma <-c(0)
  rat_new$alpha[rat_new$data_type=="traps"] <- alpha1.0
  rat_new$alpha[rat_new$data_type=="plates"] <- alpha2.0
  rat_new$alpha[rat_new$data_type=="burrows"] <- alpha3.0
  rat_new$alpha[rat_new$data_type=="faeces"] <- alpha4.0
  rat_new$alpha[rat_new$data_type=="trails"] <- alpha5.0
  
  rat_new$sigma[rat_new$data_type=="traps"] <- sigma1.0
  rat_new$sigma[rat_new$data_type=="plates"] <- sigma2.0
  rat_new$sigma[rat_new$data_type=="burrows"] <- sigma3.0
  rat_new$sigma[rat_new$data_type=="faeces"] <- sigma4.0
  rat_new$sigma[rat_new$data_type=="trails"] <- sigma5.0
  
  rat_new$lp <- rat_new$alpha+rat_new$sigma*rat_new$Rj
  rat_new$p_q_lambda <- c()
  rat_new$p_q_lamda[rat_new$data_type=="traps"] <- exp(rat_new$lp[rat_new$data_type=="traps"])
  rat_new$p_q_lamda[rat_new$data_type=="plates"] <- inv.logit(rat_new$lp[rat_new$data_type=="plates"])
  rat_new$p_q_lamda[rat_new$data_type=="burrows"] <- exp(rat_new$lp[rat_new$data_type=="burrows"])
  rat_new$p_q_lamda[rat_new$data_type=="faeces"] <- inv.logit(rat_new$lp[rat_new$data_type=="faeces"])
  rat_new$p_q_lamda[rat_new$data_type=="trails"] <- inv.logit(rat_new$lp[rat_new$data_type=="trails"])
  
  # generate data from p, q and lambda
  rat_new$t  <- 1
  rat_new$t[rat_new$offset_req==1] <- 0.5
  rat_new$outcome[rat_new$data_type=="traps"] <- rbinom(nrow(rat_new[rat_new$data_type=="traps",]),1,prob=1-exp(-rat_new$t[rat_new$data_type=="traps"]*rat_new$p_q_lamda[rat_new$data_type=="traps"]))
  rat_new$outcome[rat_new$data_type=="plates"]  <- rbinom(nrow(rat_new[rat_new$data_type=="plates",]),rat_new$offset[rat_new$data_type=="plates"], prob=rat_new$p_q_lamda[rat_new$data_type=="plates"])
  rat_new$outcome[rat_new$data_type=="burrows"] <- rpois(nrow(rat_new[rat_new$data_type=="burrows",]),lambda=rat_new$p_q_lamda[rat_new$data_type=="burrows"])
  rat_new$outcome[rat_new$data_type=="faeces"] <- rbinom(nrow(rat_new[rat_new$data_type=="faeces",]),1,prob=rat_new$p_q_lamda[rat_new$data_type=="faeces"])
  rat_new$outcome[rat_new$data_type=="trails"] <- rbinom(nrow(rat_new[rat_new$data_type=="trails",]),1,prob=rat_new$p_q_lamda[rat_new$data_type=="trails"])
  
  # Create new simulated human data
  human_new <- human
  human_new$Rj <- Rj[ID[id.human]]
  mu0.human <- as.numeric(D.human%*%gamma0)
  
  # xis
  eta.human <- mu0.human
  
  if(n.xi==1){
    human_new$lp <- eta.human + xi0*human_new$Rj
  }
  if(n.xi==2){
    human_new$lp[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*human_new$Rj[ppl.in.xi1]
    human_new$lp[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*human_new$Rj[ppl.in.xi2]
  }
  if(n.xi==3){
    human_new$lp[ppl.in.xi1] <- eta.human[ppl.in.xi1] + xi0[1]*human_new$Rj[ppl.in.xi1]
    human_new$lp[ppl.in.xi2] <- eta.human[ppl.in.xi2] + xi0[2]*human_new$Rj[ppl.in.xi2]
    human_new$lp[ppl.in.xi3] <- eta.human[ppl.in.xi3] + xi0[3]*human_new$Rj[ppl.in.xi3]
  }
  
  human_new$p <- inv.logit(human_new$lp)
  human_new$outcome <- rbinom(nrow(human_new),1,prob=human_new$p)
  
  human <- human_new
  rat <- rat_new
  
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
      
      # Human
      eta.human <- val$mu.human
      
      if(n.xi==1){
        eta.human <- eta.human + val$xi*R[ID[id.human]]
      }
      if(n.xi==2){
        eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + val$xi[1]*R[ID[id.human][ppl.in.xi1]] 
        eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + val$xi[2]*R[ID[id.human][ppl.in.xi2]]
      }
      if(n.xi==3){
        eta.human[ppl.in.xi1] <- eta.human[ppl.in.xi1] + val$xi[1]*R[ID[id.human][ppl.in.xi1]] 
        eta.human[ppl.in.xi2] <- eta.human[ppl.in.xi2] + val$xi[2]*R[ID[id.human][ppl.in.xi2]]
        eta.human[ppl.in.xi3] <- eta.human[ppl.in.xi3] + val$xi[3]*R[ID[id.human][ppl.in.xi3]]
      }
      
      llik.human <- sum(z*eta.human-log(1+exp(eta.human)))
      
      diff.R <- R-val$mu
      out <- as.numeric(-0.5*(val$log.det.Sigma+t(diff.R)%*%val$Sigma.inv%*%diff.R))+ 
        llik1+llik2+llik3+llik4+llik5+llik.human
    }
    
    n.samples <- (n.sim-burnin)/thin
    
    sim <- langevin.MCMC.samples(N, n.sim, n.samples, 
                                 estim, burnin, thin, Sigma.sroot)
    
    par0 <- name.it(with_human_S, with_human_N, with_Ui,
                    alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                    sigma1.0,sigma2.0,sigma3.0,sigma4.0,sigma5.0,
                    beta0,phi0,gamma0,xi0,omega2.0,zeta0,omega2_nugg0, psi0,
                    D.aux, glm.fit, n.xi)
    
    
    
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
      
      gamma <- par[(p+12):(p+p.human+11)]
      val$xi <- par[(p+p.human+12):(p+p.human+12+n.xi-1)]
      val$mu.human <- as.numeric(D.human%*%gamma)
      
      if(with_human_S==TRUE & with_human_N==FALSE){
        omega2 <- exp(par[p+p.human+13+n.xi-1])
        zeta <- exp(par[p+p.human+14+n.xi-1])
        ifelse(with_Ui == TRUE, psi <- exp(par[p+p.human+15+n.xi-1])/(1+exp(par[p+p.human+15+n.xi-1])), psi <- 1)
      }
      
      if(with_human_S==FALSE & with_human_N==TRUE){
        omega2_nugg <- exp(par[p+p.human+13+n.xi-1])
        ifelse(with_Ui == TRUE, psi <- exp(par[p+p.human+14+n.xi-1])/(1+exp(par[p+p.human+14+n.xi-1])), psi <- 1)
      }
      
      if(with_human_S==TRUE & with_human_N==TRUE){
        omega2 <- exp(par[p+p.human+13+n.xi-1])
        zeta <- exp(par[p+p.human+14+n.xi-1])
        omega2_nugg <- exp(par[p+p.human+15+n.xi-1])
        ifelse(with_Ui == TRUE, psi <- exp(par[p+p.human+16+n.xi-1])/(1+exp(par[p+p.human+16+n.xi-1])), psi <- 1)
      }
      
      if(with_human_S==FALSE & with_human_N==FALSE){
        ifelse(with_Ui == TRUE, psi <- exp(par[p+p.human+13+n.xi-1])/(1+exp(par[p+p.human+13+n.xi-1])), psi <- 1)
      }
      
      Sigma <- cov.matrix.setup(psi, U, phi, with_human_N, with_human_S, 
                                loc.in.xi, loc.in.human, 
                                xi0 = val$xi, omega2_nugg, omega2, zeta,
                                n.xi)
      
      val$Sigma.inv <- solve(Sigma)
      val$log.det.Sigma <- determinant(Sigma)$modulus
      
      sapply(1:(dim(sim)[1]),function(i) log.integrand(sim[i,],val))
    }
    
    log.f.tilde <- compute.log.f(par0)
    
    MC.log.lik <- function(par) {
      log(mean(exp(compute.log.f(par)-log.f.tilde)))
    }
    
    MC.log.lik(par0)
    
    estim.par <- nlminb(par0,
                        function(x) -MC.log.lik(x),
                        control=list(trace=1, rel.tol = rel.tol, iter.max = iter.max))
    
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
      message("Iteration #", k-1, " completed. Euclid norm = ", euclid.norm[k],". Total time elapsed: ", round(time.diff,3), units(time.diff), ". MCMC samples: ", n.samples)
      message("Parameter estimates: ")
      message(paste(names(par0), round(par0, 5), " | "))
    }else{
      message("Iteration #", k-1, " completed."," Total time elapsed: ", round(time.diff,3), units(time.diff), ". MCMC samples: ", n.samples)
      message("Parameter estimates: ")
      message(paste(names(par0), round(par0, 5), " | "))
    }
    
  }
  
  par.matrix <- par.matrix[, colSums(is.na(par.matrix)) != nrow(par.matrix)]
  par_rel <- par_rel[, colSums(is.na(par_rel)) != nrow(par_rel)]
  
  # unstandardised regression coeffcients
  par_us_human <- unstandardise.coeff.fn(par0[paste0("human__",names(coef(glm.fit)))],intercept = TRUE, x = D.human.us)
  par_us_rat <- unstandardise.coeff.fn(par0[paste0("rat__",c(colnames(D.aux)))],intercept = FALSE, x = D.unscale)
  par_us <- par0
  par_us[c(paste0("human__",names(coef(glm.fit))), paste0("rat__",c(colnames(D.aux))))] <- c(par_us_human,par_us_rat)
  par.list <- list(estim.par, par_us) # REDUCED TO SAVE MEMORY
  names(par.list) <- c("final estimates", "rescaled regression coefficients")
  
  return(par.list)  
}
