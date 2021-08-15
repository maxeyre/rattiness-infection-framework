# Functions for Rattiness - epidemiology project

# 1. Log-likelihood function for rattiness ecology model (without human data) with no spatial correlation
llik.non.spatial.rattiness.fn <- function(par) {
  # Set up parameters
  alpha <- par[1:n.inst]
  sigma <- exp(par[((n.inst+1):(2*n.inst))])
  n.loc <- length(unique(ID))
  
  out <- 0
  for(j in 1:n.loc) { # we use j to denote location
    
    R.halton <- U.halton
    rat.j <- rat[ID==j,] 
    data_type.j <- rat.j$data_type 
    ind.traps.j <- which(data_type.j=="traps")
    ind.plates.j <- which(data_type.j=="plates")
    ind.burrows.j <- which(data_type.j=="burrows")
    ind.faeces.j <- which(data_type.j=="faeces")
    ind.trails.j <- which(data_type.j=="trails")
    
    out.traps.j <- 1
    out.plates.j <- 1
    out.burrows.j <- 1
    out.faeces.j <- 1
    out.trails.j <- 1
    
    
    if(length(ind.traps.j) > 0) { 
      eta.traps.j <- sapply(R.halton,function(R.j) alpha[1]+sigma[1]*R.j) 
      lambda.traps.j <- exp(eta.traps.j) 
      t.j <- rep(1,length(ind.traps.j)) 
      if(any(rat.j$offset_req[ind.traps.j]==1)) t.j[rat.j$offset_req[ind.traps.j]==1] <- 0.5
      
      llik.traps.j <- sapply(1:length(ind.traps.j), 
                             function(h) dbinom(rat.j$outcome[ind.traps.j][h],1, 
                                                prob = 1-exp(-t.j[h]*lambda.traps.j),
                                                log=TRUE))
      if(length(ind.traps.j) > 1) {
        out.traps.j <- apply(llik.traps.j,1,function(row) exp(sum(row))) 
      } else {
        out.traps.j <- exp(llik.traps.j)
      }
    }
    
    if(length(ind.plates.j) > 0) {
      eta.plates.j <- sapply(R.halton,function(R.j) alpha[2]+sigma[2]*R.j)
      q.plates.j <- 1/(1+exp(-eta.plates.j)) 
      llik.plates.j <- sapply(1:length(ind.plates.j), 
                              function(h) dbinom(rat.j$outcome[ind.plates.j][h],
                                                 rat.j$offset[ind.plates.j][h],
                                                 prob = q.plates.j,
                                                 log=TRUE))
      if(length(ind.plates.j) > 1) {
        out.plates.j <- apply(llik.plates.j,1,function(row) exp(sum(row)))
      } else {
        out.plates.j <- exp(llik.plates.j)
      }
    }
    
    if(length(ind.burrows.j) > 0) { 
      eta.burrows.j <- sapply(R.halton,function(R.j) alpha[3]+sigma[3]*R.j) 
      lambda.burrows.j <- exp(eta.burrows.j) 
      
      llik.burrows.j <- sapply(1:length(ind.burrows.j), 
                               function(h) -lambda.burrows.j + rat.j$outcome[ind.burrows.j][h]*log(lambda.burrows.j))
      
      if(length(ind.burrows.j) > 1) {
        out.burrows.j <- apply(llik.burrows.j,1,function(row) exp(sum(row))) 
      } else {
        out.burrows.j <- exp(llik.burrows.j)
      }
    }
    
    if(length(ind.faeces.j)>0) {
      eta.faeces.j <- sapply(R.halton,function(R.j) alpha[4]+sigma[4]*R.j) 
      p.faeces.j <- 1/(1+exp(- eta.faeces.j)) 
      
      llik.faeces.j <- sapply(1:length(ind.faeces.j), function(h) rat.j$outcome[ind.faeces.j][h]*log(p.faeces.j)+
                                (1-rat.j$outcome[ind.faeces.j][h])*log(1-p.faeces.j))
      
      if(length(ind.faeces.j) > 1) {
        out.faeces.j <- apply(llik.faeces.j,1,function(row) exp(sum(row))) 
      } else {
        out.faeces.j <- exp(llik.faeces.j)
      }
    }
    
    if(length(ind.trails.j)>0) {
      eta.trails.j <- sapply(R.halton,function(R.j) alpha[5]+sigma[5]*R.j) 
      p.trails.j <- 1/(1+exp(- eta.trails.j)) 
      
      llik.trails.j <- sapply(1:length(ind.trails.j), function(h) rat.j$outcome[ind.trails.j][h]*log(p.trails.j)+
                                (1-rat.j$outcome[ind.trails.j][h])*log(1-p.trails.j))
      
      if(length(ind.trails.j) > 1) {
        out.trails.j <- apply(llik.trails.j,1,function(row) exp(sum(row))) 
      } else {
        out.trails.j <- exp(llik.trails.j)
      }
    }
    
    out <- out+log(mean(exp(log(out.traps.j)+log(out.plates.j)+
                              log(out.burrows.j)+log(out.faeces.j)+log(out.trails.j))))
    
  }
  return(out)
}

# 2. Compute predicted Uj values at each location (no covariates)
compute.pred.mean.Uj.no.cov.fn <- function(par,j) {
  
  alpha <- par[1:n.inst]
  sigma <- exp(par[((n.inst+1):(2*n.inst))])
  n.loc <- length(unique(ID))
  
  eta.R.j <- sum(D[j,]*beta)
  R.halton <- eta.R.j+U.halton
  
  rat.j <- rat[ID==j,]
  data_type.j <- rat.j$data_type
  ind.traps.j <- which(data_type.j=="traps")
  ind.plates.j <- which(data_type.j=="plates")
  ind.burrows.j <- which(data_type.j=="burrows")
  ind.faeces.j <- which(data_type.j=="faeces")
  ind.trails.j <- which(data_type.j=="trails")
  
  out.traps.j <- 1
  out.plates.j <- 1
  out.burrows.j <- 1
  out.faeces.j <- 1
  out.trails.j <- 1
  
  if(length(ind.traps.j) > 0) { 
    eta.traps.j <- sapply(R.halton,function(R.j) alpha[1]+sigma[1]*R.j) 
    lambda.traps.j <- exp(eta.traps.j) 
    t.j <- rep(1,length(ind.traps.j)) 
    if(any(rat.j$offset_req[ind.traps.j]==1)) t.j[rat.j$offset_req[ind.traps.j]==1] <- 0.5
    
    llik.traps.j <- sapply(1:length(ind.traps.j), 
                           function(h) dbinom(rat.j$outcome[ind.traps.j][h],1, 
                                              prob = 1-exp(-t.j[h]*lambda.traps.j),
                                              log=TRUE))
    if(length(ind.traps.j) > 1) {
      out.traps.j <- apply(llik.traps.j,1,function(row) exp(sum(row))) 
    } else {
      out.traps.j <- exp(llik.traps.j)
    }
  }
  
  if(length(ind.plates.j) > 0) {
    eta.plates.j <- sapply(R.halton,function(R.j) alpha[2]+sigma[2]*R.j)
    q.plates.j <- 1/(1+exp(-eta.plates.j)) 
    llik.plates.j <- sapply(1:length(ind.plates.j), 
                            function(h) dbinom(rat.j$outcome[ind.plates.j][h],
                                               rat.j$offset[ind.plates.j][h],
                                               prob = q.plates.j,
                                               log=TRUE))
    if(length(ind.plates.j) > 1) {
      out.plates.j <- apply(llik.plates.j,1,function(row) exp(sum(row)))
    } else {
      out.plates.j <- exp(llik.plates.j)
    }
  }
  
  if(length(ind.burrows.j) > 0) { 
    eta.burrows.j <- sapply(R.halton,function(R.j) alpha[3]+sigma[3]*R.j) 
    lambda.burrows.j <- exp(eta.burrows.j) 
    
    llik.burrows.j <- sapply(1:length(ind.burrows.j), 
                             function(h) -lambda.burrows.j + rat.j$outcome[ind.burrows.j][h]*log(lambda.burrows.j))
    
    if(length(ind.burrows.j) > 1) {
      out.burrows.j <- apply(llik.burrows.j,1,function(row) exp(sum(row))) 
    } else {
      out.burrows.j <- exp(llik.burrows.j)
    }
  }
  
  if(length(ind.faeces.j)>0) {
    eta.faeces.j <- sapply(R.halton,function(R.j) alpha[4]+sigma[4]*R.j) 
    p.faeces.j <- 1/(1+exp(- eta.faeces.j)) 
    
    llik.faeces.j <- sapply(1:length(ind.faeces.j), function(h) rat.j$outcome[ind.faeces.j][h]*log(p.faeces.j)+
                              (1-rat.j$outcome[ind.faeces.j][h])*log(1-p.faeces.j))
    
    if(length(ind.faeces.j) > 1) {
      out.faeces.j <- apply(llik.faeces.j,1,function(row) exp(sum(row))) 
    } else {
      out.faeces.j <- exp(llik.faeces.j)
    }
  }
  
  if(length(ind.trails.j)>0) {
    eta.trails.j <- sapply(R.halton,function(R.j) alpha[5]+sigma[5]*R.j) 
    p.trails.j <- 1/(1+exp(- eta.trails.j)) 
    
    llik.trails.j <- sapply(1:length(ind.trails.j), function(h) rat.j$outcome[ind.trails.j][h]*log(p.trails.j)+
                              (1-rat.j$outcome[ind.trails.j][h])*log(1-p.trails.j))
    
    if(length(ind.trails.j) > 1) {
      out.trails.j <- apply(llik.trails.j,1,function(row) exp(sum(row))) 
    } else {
      out.trails.j <- exp(llik.trails.j)
    }
  }
  
  den <- mean(exp(log(out.traps.j)+log(out.plates.j)+
                    log(out.burrows.j)+log(out.faeces.j)+log(out.trails.j)))
  num <- mean(U.halton*exp(log(out.traps.j)+log(out.plates.j)+
                             log(out.burrows.j)+log(out.faeces.j)+log(out.trails.j)))
  num/den
}

# 3. Fit non-spatial rattiness-ecology model with covariates
llik.non.spatial.rattiness.with.covariates.fn <- function(par) {
  # Set up parameters
  alpha <- par[1:n.inst]
  sigma <- exp(par[((n.inst+1):(2*n.inst))])
  n.loc <- length(unique(ID))
  
  beta <- par[(2*n.inst+1):(2*n.inst+p)]
  eta.R <- as.numeric(D%*%beta)
  
  out <- 0
  for(j in 1:n.loc) { # we use j to denote location
    
    R.halton <- eta.R[j] + U.halton
    rat.j <- rat[ID==j,] 
    data_type.j <- rat.j$data_type 
    ind.traps.j <- which(data_type.j=="traps")
    ind.plates.j <- which(data_type.j=="plates")
    ind.burrows.j <- which(data_type.j=="burrows")
    ind.faeces.j <- which(data_type.j=="faeces")
    ind.trails.j <- which(data_type.j=="trails")
    
    out.traps.j <- 1
    out.plates.j <- 1
    out.burrows.j <- 1
    out.faeces.j <- 1
    out.trails.j <- 1
    
    
    if(length(ind.traps.j) > 0) { 
      eta.traps.j <- sapply(R.halton,function(R.j) alpha[1]+sigma[1]*R.j) 
      lambda.traps.j <- exp(eta.traps.j) 
      t.j <- rep(1,length(ind.traps.j)) 
      if(any(rat.j$offset_req[ind.traps.j]==1)) t.j[rat.j$offset_req[ind.traps.j]==1] <- 0.5
      
      llik.traps.j <- sapply(1:length(ind.traps.j), 
                             function(h) dbinom(rat.j$outcome[ind.traps.j][h],1, 
                                                prob = 1-exp(-t.j[h]*lambda.traps.j),
                                                log=TRUE))
      if(length(ind.traps.j) > 1) {
        out.traps.j <- apply(llik.traps.j,1,function(row) exp(sum(row))) 
      } else {
        out.traps.j <- exp(llik.traps.j)
      }
    }
    
    if(length(ind.plates.j) > 0) {
      eta.plates.j <- sapply(R.halton,function(R.j) alpha[2]+sigma[2]*R.j)
      q.plates.j <- 1/(1+exp(-eta.plates.j)) 
      llik.plates.j <- sapply(1:length(ind.plates.j), 
                              function(h) dbinom(rat.j$outcome[ind.plates.j][h],
                                                 rat.j$offset[ind.plates.j][h],
                                                 prob = q.plates.j,
                                                 log=TRUE))
      if(length(ind.plates.j) > 1) {
        out.plates.j <- apply(llik.plates.j,1,function(row) exp(sum(row)))
      } else {
        out.plates.j <- exp(llik.plates.j)
      }
    }
    
    if(length(ind.burrows.j) > 0) { 
      eta.burrows.j <- sapply(R.halton,function(R.j) alpha[3]+sigma[3]*R.j) 
      lambda.burrows.j <- exp(eta.burrows.j) 
      
      llik.burrows.j <- sapply(1:length(ind.burrows.j), 
                               function(h) -lambda.burrows.j + rat.j$outcome[ind.burrows.j][h]*log(lambda.burrows.j))
      
      if(length(ind.burrows.j) > 1) {
        out.burrows.j <- apply(llik.burrows.j,1,function(row) exp(sum(row))) 
      } else {
        out.burrows.j <- exp(llik.burrows.j)
      }
    }
    
    if(length(ind.faeces.j)>0) {
      eta.faeces.j <- sapply(R.halton,function(R.j) alpha[4]+sigma[4]*R.j) 
      p.faeces.j <- 1/(1+exp(- eta.faeces.j)) 
      
      llik.faeces.j <- sapply(1:length(ind.faeces.j), function(h) rat.j$outcome[ind.faeces.j][h]*log(p.faeces.j)+
                                (1-rat.j$outcome[ind.faeces.j][h])*log(1-p.faeces.j))
      
      if(length(ind.faeces.j) > 1) {
        out.faeces.j <- apply(llik.faeces.j,1,function(row) exp(sum(row))) 
      } else {
        out.faeces.j <- exp(llik.faeces.j)
      }
    }
    
    if(length(ind.trails.j)>0) {
      eta.trails.j <- sapply(R.halton,function(R.j) alpha[5]+sigma[5]*R.j) 
      p.trails.j <- 1/(1+exp(- eta.trails.j)) 
      
      llik.trails.j <- sapply(1:length(ind.trails.j), function(h) rat.j$outcome[ind.trails.j][h]*log(p.trails.j)+
                                (1-rat.j$outcome[ind.trails.j][h])*log(1-p.trails.j))
      
      if(length(ind.trails.j) > 1) {
        out.trails.j <- apply(llik.trails.j,1,function(row) exp(sum(row))) 
      } else {
        out.trails.j <- exp(llik.trails.j)
      }
    }
    
    out <- out+log(mean(exp(log(out.traps.j)+log(out.plates.j)+
                              log(out.burrows.j)+log(out.faeces.j)+log(out.trails.j))))
    
  }
  return(out)
}

# 4. Compute predicted Uj values at each location (with covariates)
# compute.pred.mean.Uj.with.cov.fn <- function(par,j) {
#   n.loc <- length(unique(ID))
#   
#   alpha <- par[1:n.inst]
#   sigma <- exp(par[((n.inst+1):(2*n.inst))])
#   beta <- par[(2*n.inst+1):(2*n.inst+p)]
#   
#   n.loc <- length(unique(ID))
#   
#   eta.R.j <- sum(D[j,]*beta)
#   
#   R.halton <- eta.R.j+U.halton
#   R.halton <- U.halton
#   rat.j <- rat[ID==j,]
#   data_type.j <- rat.j$data_type
#   ind.plates.j <- which(data_type.j=="plates")
#   ind.traps.j <- which(data_type.j=="traps")
#   ind.signs.j <- which(data_type.j=="signs")
#   
#   out.signs.j <- 1
#   out.traps.j <- 1
#   out.plates.j <- 1
#   
#   if(length(ind.signs.j)>0) {
#     eta.signs.j <- sapply(R.halton,function(R.j) alpha[1]+sigma[1]*R.j)
#     p.signs.j <- 1/(1+exp(- eta.signs.j))
#     
#     llik.signs.j <- 
#       rat.j$outcome[ind.signs.j]*log(p.signs.j)+
#       (1-rat.j$outcome[ind.signs.j])*log(1-p.signs.j)
#     out.signs.j <- exp(llik.signs.j)
#   }
#   
#   if(length(ind.traps.j) > 0) {
#     eta.traps.j <- sapply(R.halton,function(R.j) alpha[2]+sigma[2]*R.j)
#     lambda.traps.j <- exp(eta.traps.j)
#     t.j <- rep(1,length(ind.traps.j))
#     if(any(rat.j$offset_req[ind.traps.j]==1)) t.j[rat.j$offset_req[ind.traps.j]==1] <- 0.5
#     llik.traps.j <- sapply(1:length(ind.traps.j), 
#                            function(h) dbinom(rat.j$outcome[ind.traps.j][h],1,
#                                               prob = 1-exp(-t.j[h]*lambda.traps.j),
#                                               log=TRUE))
#     if(length(ind.traps.j) > 1) {
#       out.traps.j <- apply(llik.traps.j,1,function(row) exp(sum(row)))
#     } else {
#       out.traps.j <- exp(llik.traps.j)
#     }
#   }
#   
#   if(length(ind.plates.j) > 0) {
#     eta.plates.j <- sapply(R.halton,function(R.j) alpha[3]+sigma[3]*R.j)
#     q.plates.j <- 1/(1+exp(-eta.plates.j))
#     llik.plates.j <- sapply(1:length(ind.plates.j), 
#                             function(h) dbinom(rat.j$outcome[ind.plates.j][h],
#                                                rat.j$offset[ind.plates.j][h],
#                                                prob = q.plates.j,
#                                                log=TRUE))
#     if(length(ind.plates.j) > 1) {
#       out.plates.j <- apply(llik.plates.j,1,function(row) exp(sum(row)))
#     } else {
#       out.plates.j <- exp(llik.plates.j)
#     }
#   }
#   den <- mean(exp(log(out.signs.j)+log(out.traps.j)+log(out.plates.j)))
#   num <- mean(U.halton*exp(log(out.signs.j)+log(out.traps.j)+log(out.plates.j)))
#   num/den
# }

# compute.pred.mean.Uj.with.cov.fn <- function(par,j) {
#   n.loc <- length(unique(ID))
#   
#   alpha <- par[1:n.inst]
#   sigma <- exp(par[((n.inst+1):(2*n.inst))])
#   beta <- par[(2*n.inst+1):(2*n.inst+p)]
#   
#   n.loc <- length(unique(ID))
#   
#   eta.R.j <- sum(D[j,]*beta)
#   
#   R.halton <- eta.R.j+U.halton
#   R.halton <- U.halton
#   rat.j <- rat[ID==j,]
#   data_type.j <- rat.j$data_type
#   ind.traps.j <- which(data_type.j=="traps")
#   ind.plates.j <- which(data_type.j=="plates")
#   ind.burrows.j <- which(data_type.j=="burrows")
#   ind.faeces.j <- which(data_type.j=="faeces")
#   ind.trails.j <- which(data_type.j=="trails")
#   
#   out.traps.j <- 1
#   out.plates.j <- 1
#   out.burrows.j <- 1
#   out.faeces.j <- 1
#   out.trails.j <- 1
#   
#   if(length(ind.traps.j) > 0) { 
#     eta.traps.j <- sapply(R.halton,function(R.j) alpha[1]+sigma[1]*R.j) 
#     lambda.traps.j <- exp(eta.traps.j) 
#     t.j <- rep(1,length(ind.traps.j)) 
#     if(any(rat.j$offset_req[ind.traps.j]==1)) t.j[rat.j$offset_req[ind.traps.j]==1] <- 0.5
#     
#     llik.traps.j <- sapply(1:length(ind.traps.j), 
#                            function(h) dbinom(rat.j$outcome[ind.traps.j][h],1, 
#                                               prob = 1-exp(-t.j[h]*lambda.traps.j),
#                                               log=TRUE))
#     if(length(ind.traps.j) > 1) {
#       out.traps.j <- apply(llik.traps.j,1,function(row) exp(sum(row))) 
#     } else {
#       out.traps.j <- exp(llik.traps.j)
#     }
#   }
#   
#   if(length(ind.plates.j) > 0) {
#     eta.plates.j <- sapply(R.halton,function(R.j) alpha[2]+sigma[2]*R.j)
#     q.plates.j <- 1/(1+exp(-eta.plates.j)) 
#     llik.plates.j <- sapply(1:length(ind.plates.j), 
#                             function(h) dbinom(rat.j$outcome[ind.plates.j][h],
#                                                rat.j$offset[ind.plates.j][h],
#                                                prob = q.plates.j,
#                                                log=TRUE))
#     if(length(ind.plates.j) > 1) {
#       out.plates.j <- apply(llik.plates.j,1,function(row) exp(sum(row)))
#     } else {
#       out.plates.j <- exp(llik.plates.j)
#     }
#   }
#   
#   if(length(ind.burrows.j) > 0) { 
#     eta.burrows.j <- sapply(R.halton,function(R.j) alpha[3]+sigma[3]*R.j) 
#     lambda.burrows.j <- exp(eta.burrows.j) 
#     
#     llik.burrows.j <- sapply(1:length(ind.burrows.j), 
#                              function(h) -lambda.burrows.j + rat.j$outcome[ind.burrows.j][h]*log(lambda.burrows.j))
#     
#     if(length(ind.burrows.j) > 1) {
#       out.burrows.j <- apply(llik.burrows.j,1,function(row) exp(sum(row))) 
#     } else {
#       out.burrows.j <- exp(llik.burrows.j)
#     }
#   }
#   
#   if(length(ind.faeces.j)>0) {
#     eta.faeces.j <- sapply(R.halton,function(R.j) alpha[4]+sigma[4]*R.j) 
#     p.faeces.j <- 1/(1+exp(- eta.faeces.j)) 
#     
#     llik.faeces.j <- sapply(1:length(ind.faeces.j), function(h) rat.j$outcome[ind.faeces.j][h]*log(p.faeces.j)+
#                               (1-rat.j$outcome[ind.faeces.j][h])*log(1-p.faeces.j))
#     
#     if(length(ind.faeces.j) > 1) {
#       out.faeces.j <- apply(llik.faeces.j,1,function(row) exp(sum(row))) 
#     } else {
#       out.faeces.j <- exp(llik.faeces.j)
#     }
#   }
#   
#   if(length(ind.trails.j)>0) {
#     eta.trails.j <- sapply(R.halton,function(R.j) alpha[5]+sigma[5]*R.j) 
#     p.trails.j <- 1/(1+exp(- eta.trails.j)) 
#     
#     llik.trails.j <- sapply(1:length(ind.trails.j), function(h) rat.j$outcome[ind.trails.j][h]*log(p.trails.j)+
#                               (1-rat.j$outcome[ind.trails.j][h])*log(1-p.trails.j))
#     
#     if(length(ind.trails.j) > 1) {
#       out.trails.j <- apply(llik.trails.j,1,function(row) exp(sum(row))) 
#     } else {
#       out.trails.j <- exp(llik.trails.j)
#     }
#   }
#   
#   den <- mean(exp(log(out.traps.j)+log(out.plates.j)+
#                     log(out.burrows.j)+log(out.faeces.j)+log(out.trails.j)))
#   num <- mean(U.halton*exp(log(out.traps.j)+log(out.plates.j)+
#                              log(out.burrows.j)+log(out.faeces.j)+log(out.trails.j)))
#   num/den
# }
