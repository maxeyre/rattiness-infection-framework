cov.matrix.setup <- function(psi0, U, phi0, with_human_N, with_human_S, 
                             loc.in.xi, loc.in.human, 
                             xi0, omega2_nugg0, omega2.0, zeta0,
                             n.xi){
  Sigma0 <- as.matrix(psi0*exp(-U/phi0)) # cov(rat,rat)
  diag(Sigma0) <- 1
  
  if(n.xi==1){
    if(with_human_S==TRUE){
      Sigma0[loc.in.human,loc.in.human] <- Sigma0[loc.in.human,loc.in.human]+(1/(xi0^2))*omega2.0*exp(-U[loc.in.human,loc.in.human]/zeta0)}
    if(with_human_N==TRUE){
      diag(Sigma0[loc.in.human,loc.in.human]) <- diag(Sigma0[loc.in.human,loc.in.human]) + (1/(xi0^2))*omega2_nugg0}
  }
  
  if(n.xi==2){
    if(with_human_S==TRUE){ # (with human Gaussian process)
      
      # human_low and human_low
      Sigma0[loc.in.xi[[1]],loc.in.xi[[1]]] <- Sigma0[loc.in.xi[[1]],loc.in.xi[[1]]] + (1/(xi0[1]^2))*omega2.0*exp(-U[loc.in.xi[[1]],loc.in.xi[[1]]]/zeta0)
      # human_high and human_high
      Sigma0[loc.in.xi[[2]],loc.in.xi[[2]]] <- Sigma0[loc.in.xi[[2]],loc.in.xi[[2]]] + (1/(xi0[2]^2))*omega2.0*exp(-U[loc.in.xi[[2]],loc.in.xi[[2]]]/zeta0)
      # human_low and human_high
      Sigma0[loc.in.xi[[1]],loc.in.xi[[2]]] <- Sigma0[loc.in.xi[[1]],loc.in.xi[[2]]] + (1/(xi0[1]*xi0[2]))*omega2.0*exp(-U[loc.in.xi[[1]],loc.in.xi[[2]]]/zeta0)
      # human_high and human_low
      Sigma0[loc.in.xi[[2]],loc.in.xi[[1]]] <- Sigma0[loc.in.xi[[2]],loc.in.xi[[1]]] + (1/(xi0[2]*xi0[1]))*omega2.0*exp(-U[loc.in.xi[[2]],loc.in.xi[[1]]]/zeta0)
    }
    
    if(with_human_N==TRUE){ # (with human nugget effect)
      # human_low and human_low
      diag(Sigma0[loc.in.xi[[1]],loc.in.xi[[1]]]) <- diag(Sigma0[loc.in.xi[[1]],loc.in.xi[[1]]]) + (1/(xi0[1]^2))*omega2_nugg0
      # human_high and human_high
      diag(Sigma0[loc.in.xi[[2]],loc.in.xi[[2]]]) <- diag(Sigma0[loc.in.xi[[2]],loc.in.xi[[2]]]) + (1/(xi0[2]^2))*omega2_nugg0
    }
  }
  
  if(n.xi==3){
    if(with_human_S==TRUE){ # (with human Gaussian process)
      
      # human_1 and human_1
      Sigma0[loc.in.xi[[1]],loc.in.xi[[1]]] <- Sigma0[loc.in.xi[[1]],loc.in.xi[[1]]] + (1/(xi0[1]^2))*omega2.0*exp(-U[loc.in.xi[[1]],loc.in.xi[[1]]]/zeta0)
      # human_2 and human_2
      Sigma0[loc.in.xi[[2]],loc.in.xi[[2]]] <- Sigma0[loc.in.xi[[2]],loc.in.xi[[2]]] + (1/(xi0[2]^2))*omega2.0*exp(-U[loc.in.xi[[2]],loc.in.xi[[2]]]/zeta0)
      # human_3 and human_3
      Sigma0[loc.in.xi[[3]],loc.in.xi[[3]]] <- Sigma0[loc.in.xi[[3]],loc.in.xi[[3]]] + (1/(xi0[3]^2))*omega2.0*exp(-U[loc.in.xi[[3]],loc.in.xi[[3]]]/zeta0)
      # human_1 and human_2
      Sigma0[loc.in.xi[[1]],loc.in.xi[[2]]] <- Sigma0[loc.in.xi[[1]],loc.in.xi[[2]]] + (1/(xi0[1]*xi0[2]))*omega2.0*exp(-U[loc.in.xi[[1]],loc.in.xi[[2]]]/zeta0)
      # human_2 and human_1
      Sigma0[loc.in.xi[[2]],loc.in.xi[[1]]] <- Sigma0[loc.in.xi[[2]],loc.in.xi[[1]]] + (1/(xi0[2]*xi0[1]))*omega2.0*exp(-U[loc.in.xi[[2]],loc.in.xi[[1]]]/zeta0)
      # human_1 and human_3
      Sigma0[loc.in.xi[[1]],loc.in.xi[[3]]] <- Sigma0[loc.in.xi[[1]],loc.in.xi[[3]]] + (1/(xi0[1]*xi0[3]))*omega2.0*exp(-U[loc.in.xi[[1]],loc.in.xi[[3]]]/zeta0)
      # human_3 and human_1
      Sigma0[loc.in.xi[[3]],loc.in.xi[[1]]] <- Sigma0[loc.in.xi[[3]],loc.in.xi[[1]]] + (1/(xi0[3]*xi0[1]))*omega2.0*exp(-U[loc.in.xi[[3]],loc.in.xi[[1]]]/zeta0)
      # human_2 and human_3
      Sigma0[loc.in.xi[[2]],loc.in.xi[[3]]] <- Sigma0[loc.in.xi[[2]],loc.in.xi[[3]]] + (1/(xi0[2]*xi0[3]))*omega2.0*exp(-U[loc.in.xi[[2]],loc.in.xi[[3]]]/zeta0)
      # human_3 and human_2
      Sigma0[loc.in.xi[[3]],loc.in.xi[[2]]] <- Sigma0[loc.in.xi[[3]],loc.in.xi[[2]]] + (1/(xi0[3]*xi0[2]))*omega2.0*exp(-U[loc.in.xi[[3]],loc.in.xi[[2]]]/zeta0)
    }
    
    if(with_human_N==TRUE){ # (with human nugget effect)
      # human_1 and human_1
      diag(Sigma0[loc.in.xi[[1]],loc.in.xi[[1]]]) <- diag(Sigma0[loc.in.xi[[1]],loc.in.xi[[1]]]) + (1/(xi0[1]^2))*omega2_nugg0
      # human_2 and human_2
      diag(Sigma0[loc.in.xi[[2]],loc.in.xi[[2]]]) <- diag(Sigma0[loc.in.xi[[2]],loc.in.xi[[2]]]) + (1/(xi0[2]^2))*omega2_nugg0
      # human_3 and human_3
      diag(Sigma0[loc.in.xi[[3]],loc.in.xi[[3]]]) <- diag(Sigma0[loc.in.xi[[3]],loc.in.xi[[3]]]) + (1/(xi0[3]^2))*omega2_nugg0
    }
  }
  
  return(Sigma0)
}

name.it <- function(with_human_S, with_human_N, with_Ui,
                    alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                    sigma1.0,sigma2.0,sigma3.0,sigma4.0,sigma5.0,
                    beta0,phi0,gamma0,xi0,omega2.0,zeta0,omega2_nugg0, psi0,
                    D.aux, glm.fit, n.xi){
  if(with_human_S == TRUE & with_human_N == FALSE){
    if(with_Ui == TRUE){
      par0 <- c(alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                log(sigma1.0),log(sigma2.0),log(sigma3.0),log(sigma4.0),log(sigma5.0),
                beta0,log(phi0),gamma0,xi0,log(omega2.0),log(zeta0),log(psi0/(1-psi0)))
      names(par0) <- c("alpha_traps","alpha_plates","alpha_burrows","alpha_faeces","alpha_trails",
                       "log(sigma_traps)","log(sigma_plates)","log(sigma_burrows)","log(sigma_faeces)","log(sigma_trails)",
                       paste0("rat__",c(colnames(D.aux))),"log(phi)",paste0("human__",names(coef(glm.fit))),paste0("xi_",1:n.xi),"log(omega2)","log(zeta)","log(psi/(1-psi))")
    }else{
      par0 <- c(alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                log(sigma1.0),log(sigma2.0),log(sigma3.0),log(sigma4.0),log(sigma5.0),
                beta0,log(phi0),gamma0,xi0,log(omega2.0),log(zeta0))
      names(par0) <- c("alpha_traps","alpha_plates","alpha_burrows","alpha_faeces","alpha_trails",
                       "log(sigma_traps)","log(sigma_plates)","log(sigma_burrows)","log(sigma_faeces)","log(sigma_trails)",
                       paste0("rat__",c(colnames(D.aux))),"log(phi)",paste0("human__",names(coef(glm.fit))),paste0("xi_",1:n.xi),"log(omega2)","log(zeta)")
    }
  }
  
  if(with_human_S == FALSE & with_human_N == TRUE){
    if(with_Ui == TRUE){
      par0 <- c(alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                log(sigma1.0),log(sigma2.0),log(sigma3.0),log(sigma4.0),log(sigma5.0),
                beta0,log(phi0),gamma0,xi0,log(omega2_nugg0),log(psi0/(1-psi0)))
      names(par0) <- c("alpha_traps","alpha_plates","alpha_burrows","alpha_faeces","alpha_trails",
                       "log(sigma_traps)","log(sigma_plates)","log(sigma_burrows)","log(sigma_faeces)","log(sigma_trails)",
                       paste0("rat__",c(colnames(D.aux))),"log(phi)",paste0("human__",names(coef(glm.fit))),paste0("xi_",1:n.xi),"log(omega2_nugg)","log(psi/(1-psi))")
    }else{
      par0 <- c(alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                log(sigma1.0),log(sigma2.0),log(sigma3.0),log(sigma4.0),log(sigma5.0),
                beta0,log(phi0),gamma0,xi0,log(omega2_nugg0))
      names(par0) <- c("alpha_traps","alpha_plates","alpha_burrows","alpha_faeces","alpha_trails",
                       "log(sigma_traps)","log(sigma_plates)","log(sigma_burrows)","log(sigma_faeces)","log(sigma_trails)",
                       paste0("rat__",c(colnames(D.aux))),"log(phi)",paste0("human__",names(coef(glm.fit))),paste0("xi_",1:n.xi),"log(omega2_nugg)")
    }
  }
  
  if(with_human_S == TRUE & with_human_N == TRUE){
    if(with_Ui == TRUE){
      par0 <- c(alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                log(sigma1.0),log(sigma2.0),log(sigma3.0),log(sigma4.0),log(sigma5.0),
                beta0,log(phi0),gamma0,xi0,log(omega2.0),log(zeta0),log(omega2_nugg0),log(psi0/(1-psi0)))
      names(par0) <- c("alpha_traps","alpha_plates","alpha_burrows","alpha_faeces","alpha_trails",
                       "log(sigma_traps)","log(sigma_plates)","log(sigma_burrows)","log(sigma_faeces)","log(sigma_trails)",
                       paste0("rat__",c(colnames(D.aux))),"log(phi)",paste0("human__",names(coef(glm.fit))),paste0("xi_",1:n.xi),"log(omega2)","log(zeta)","log(omega2_nugg)","log(psi/(1-psi))")
    }else{
      par0 <- c(alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                log(sigma1.0),log(sigma2.0),log(sigma3.0),log(sigma4.0),log(sigma5.0),
                beta0,log(phi0),gamma0,xi0,log(omega2.0),log(zeta0),log(omega2_nugg0))
      names(par0) <- c("alpha_traps","alpha_plates","alpha_burrows","alpha_faeces","alpha_trails",
                       "log(sigma_traps)","log(sigma_plates)","log(sigma_burrows)","log(sigma_faeces)","log(sigma_trails)",
                       paste0("rat__",c(colnames(D.aux))),"log(phi)",paste0("human__",names(coef(glm.fit))),paste0("xi_",1:n.xi),"log(omega2)","log(zeta)","log(omega2_nugg)")
    }
  }
  
  if(with_human_S == FALSE & with_human_N == FALSE){
    if(with_Ui == TRUE){
      par0 <- c(alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                log(sigma1.0),log(sigma2.0),log(sigma3.0),log(sigma4.0),log(sigma5.0),
                beta0,log(phi0),gamma0,xi0,log(psi0/(1-psi0)))
      names(par0) <- c("alpha_traps","alpha_plates","alpha_burrows","alpha_faeces","alpha_trails",
                       "log(sigma_traps)","log(sigma_plates)","log(sigma_burrows)","log(sigma_faeces)","log(sigma_trails)",
                       paste0("rat__",c(colnames(D.aux))),"log(phi)",paste0("human__",names(coef(glm.fit))),paste0("xi_",1:n.xi),"log(psi/(1-psi))")
    }else{
      par0 <- c(alpha1.0,alpha2.0,alpha3.0,alpha4.0,alpha5.0,
                log(sigma1.0),log(sigma2.0),log(sigma3.0),log(sigma4.0),log(sigma5.0),
                beta0,log(phi0),gamma0,xi0)
      names(par0) <- c("alpha_traps","alpha_plates","alpha_burrows","alpha_faeces","alpha_trails",
                       "log(sigma_traps)","log(sigma_plates)","log(sigma_burrows)","log(sigma_faeces)","log(sigma_trails)",
                       paste0("rat__",c(colnames(D.aux))),"log(phi)",paste0("human__",names(coef(glm.fit))),paste0("xi_",1:n.xi))
    }
  }
  return(par0)
}
