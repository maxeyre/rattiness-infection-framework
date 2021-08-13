## Includes functions for standardising/unstandardising regression coefficients

# Please note that this works for a model.matrix set-up to ensure all covariates are correctly scaled.

## standardise.coeff.fn
# beta: set of regression coefficients
# intercept: logical for whether including intercept or not (TRUE = included/FALSE = not included). 
# (please note that intercept must be included if you have discrete covariates)
# x: model matrix of unscaled covariate values (unscaled) with columns in the same order as regression coefficients
standardise.coeff.fn <- function(beta, intercept=TRUE, x){
  if(is.logical(intercept)==FALSE){stop("intercept must be logical: TRUE (with intercept) or FALSE (without intercept) ")}
  
  n <- length(beta)
  
  if(intercept==TRUE){
    x.s <- scale(x[,2:n])
    center <- attr(x.s,"scaled:center")
    scaled <- attr(x.s,"scaled:scale")
    
    beta0.s <- beta[1] + sum(beta[2:n]*center)
    beta.s <- beta[2:n]*scaled
    beta.full.s <- c(beta0.s,beta.s)
  }
  
  if(intercept==FALSE){
    x.s <- scale(x)
    center <- attr(x.s,"scaled:center")
    scaled <- attr(x.s,"scaled:scale")
    
    beta.full.s <- beta[1:n]*scaled
  }
  
  return(beta.full.s)
}

## unstandardise.coeff.fn
# beta.s: set of regression coefficients estimated on standardised data
# intercept: logical for whether including intercept or not (TRUE = included/FALSE = not included). 
# (please note that intercept must be included if you have discrete covariates)
# x: model matrix of unscaled covariate values (unscaled) with columns in the same order as regression coefficients
unstandardise.coeff.fn <- function(beta.s, intercept = TRUE, x){
  if(is.logical(intercept)==FALSE){stop("intercept must be logical: TRUE (with intercept) or FALSE (without intercept) ")}
  
  n <- length(beta.s)
  
  if(intercept==TRUE){
    x.s <- scale(x[,2:n])
    center <- attr(x.s,"scaled:center")
    scaled <- attr(x.s,"scaled:scale")
    
    beta0.us <- beta.s[1] - sum(beta.s[2:n]*center/scaled)
    beta.us <- beta.s[2:n]/scaled
    beta.full.us <- c(beta0.us,beta.us)
  }
  
  if(intercept==FALSE){
    x.s <- scale(x)
    center <- attr(x.s,"scaled:center")
    scaled <- attr(x.s,"scaled:scale")
    
    beta.full.us <- beta.s[1:n]/scaled
  }
  
  return(beta.full.us)
}


# clean up of all parameter estimates and 95%CIs to unstandardise, put on log-odds scale etc.
clean.all.par.estimates.fn <- function(control, rat, human, boot_cis, scale_divide_10, scale_multiply_10, splines){
  
  sp <- function(x,k) max(0,x-k)
  sp <- Vectorize(sp)
  
  # read in controls
  multi.xi.on <- control$multi.xi.on[1]
  xi.var <- control$xi.var[1]
  cov.rat <- control$rat[!is.na(control$rat)]
  cov.rat.sub <- cov.rat[str_detect(cov.rat, "sp.", negate=TRUE) & str_detect(cov.rat, "as.factor", negate=TRUE)]
  cov.rat.sub <- unique(c(cov.rat.sub, unlist(str_extract_all(cov.rat[str_detect(cov.rat, ",", negate=TRUE)],  "(?<=\\().+?(?=\\))"))))
  cov.human <- control$human[!is.na(control$human)]
  if(multi.xi.on == TRUE){cov.human <- c(cov.human, paste0("as.factor(",xi.var,")"))}
  cov.human.sub <- cov.human[str_detect(cov.human, "sp.", negate=TRUE) & str_detect(cov.human, "as.factor", negate=TRUE)]
  cov.human.sub <- unique(c(cov.human.sub, unlist(str_extract_all(cov.human[str_detect(cov.human, ",", negate=TRUE)],  "(?<=\\().+?(?=\\))")), cov.rat.sub))
  formula.rat <- as.formula(paste("~-1+",paste(cov.rat,collapse="+")))
  formula.human <- as.formula(paste("outcome ~",paste(cov.human,collapse="+")))
  
  rat <- na.omit(rat[,c("X","Y","data_type","outcome","offset","offset_req", cov.rat.sub)])
  human <- na.omit(human[,c("X","Y","outcome",cov.human.sub)])
  
  # rattiness covariate model matrix at all locations
  rat.human <- bind_rows(rat[,cov.rat.sub], human[,cov.rat.sub])
  D.aux <- as.matrix(model.matrix(formula.rat, data=rat.human))
  
  # coordinates
  coords.rat <- as.matrix(rat[,c("X","Y")])
  coords.human <- as.matrix(human[,c("X","Y")])
  coords.set <- data.frame(rbind(coords.rat,coords.human))
  colnames(coords.set)  <- c("X","Y")
  id.rat <- 1:nrow(coords.rat)
  id.human <- (nrow(coords.rat)+1):(nrow(coords.rat)+nrow(coords.human))
  ID <- create.ID.coords(coords.set,~X+Y)
  coords <- unique(coords.set)
  
  # standardise rattiness covariates
  p <- ncol(D.aux)
  N <- nrow(coords)
  D <- matrix(NA,nrow=N,ncol=p)
  for(i in 1:p) {
    D[,i] <- tapply(D.aux[,i],ID,max)
  }
  D.unscale <- D
  
  # rat
  boot_cis$estimate[str_detect(boot_cis$parameters, "rat__")] <- unstandardise.coeff.fn(boot_cis$estimate[str_detect(boot_cis$parameters, "rat__")],intercept = FALSE, x = D.unscale)
  boot_cis$lci[str_detect(boot_cis$parameters, "rat__")] <- unstandardise.coeff.fn(boot_cis$lci[str_detect(boot_cis$parameters, "rat__")],intercept = FALSE, x = D.unscale)
  boot_cis$uci[str_detect(boot_cis$parameters, "rat__")] <- unstandardise.coeff.fn(boot_cis$uci[str_detect(boot_cis$parameters, "rat__")],intercept = FALSE, x = D.unscale)
  
  # human
  glm.fit <- glm(formula.human, data=human,family=binomial,x=TRUE)
  boot_cis$estimate[str_detect(boot_cis$parameters, "human__")] <- unstandardise.coeff.fn(boot_cis$estimate[str_detect(boot_cis$parameters, "human__")],intercept = TRUE, x = glm.fit$x)
  boot_cis$lci[str_detect(boot_cis$parameters, "human__")] <- unstandardise.coeff.fn(boot_cis$lci[str_detect(boot_cis$parameters, "human__")],intercept = TRUE, x = glm.fit$x)
  boot_cis$uci[str_detect(boot_cis$parameters, "human__")] <- unstandardise.coeff.fn(boot_cis$uci[str_detect(boot_cis$parameters, "human__")],intercept = TRUE, x = glm.fit$x)
  
  # scale_divide_10
  if(!is.null(scale_divide_10)){
    for(k in 1:length(scale_divide_10)){
      boot_cis[str_detect(boot_cis$parameters, scale_divide_10[k]),2:4] <- boot_cis[str_detect(boot_cis$parameters, scale_divide_10[k]),2:4]/10
    }
  }
  
  # scale_multiply_10
  if(!is.null(scale_multiply_10)){
    for(k in 1:length(scale_multiply_10)){
      boot_cis[str_detect(boot_cis$parameters, scale_multiply_10[k]),2:4] <- boot_cis[str_detect(boot_cis$parameters, scale_multiply_10[k]),2:4]*10
    }
  }
  
  # splines
  if(!is.null(splines)){
    for(k in 1:length(splines)){
      nk <- nrow(boot_cis[str_detect(boot_cis$parameters, splines[k]),2:4])
      boot_cis[str_detect(boot_cis$parameters, splines[k]),2:4] <- boot_cis[str_detect(boot_cis$parameters, splines[k]),2:4] + 
        cumsum(c(0,lag(boot_cis$estimate[str_detect(boot_cis$parameters, splines[k])])[2:nk]))
      }
  }
  
  # exponentiate for odds-scale
  boot_cis$estimate[str_detect(boot_cis$parameters, "human__")] <- exp(boot_cis$estimate[str_detect(boot_cis$parameters, "human__")])
  boot_cis$lci[str_detect(boot_cis$parameters, "human__")] <- exp(boot_cis$lci[str_detect(boot_cis$parameters, "human__")])
  boot_cis$uci[str_detect(boot_cis$parameters, "human__")] <- exp(boot_cis$uci[str_detect(boot_cis$parameters, "human__")])
  
  # sigma parameters
  boot_cis[str_detect(boot_cis$parameters,"log\\(sigma_"),2:4] <- exp(boot_cis[str_detect(boot_cis$parameters,"log\\(sigma_"),2:4])
  boot_cis[str_detect(boot_cis$parameters,"log\\(sigma_"),1] <- sapply(boot_cis[str_detect(boot_cis$parameters,"log\\(sigma_"),1], function(x) substr(x,5, nchar(x)-1))
  
  # phi, psi, omega2
  boot_cis[boot_cis$parameters=="log(phi)",2:4] <- exp(boot_cis[boot_cis$parameters=="log(phi)",2:4])
  boot_cis[boot_cis$parameters=="log(phi)",1] <- "phi"
  boot_cis[boot_cis$parameters=="log(psi/(1-psi))",2:4] <- exp(boot_cis[boot_cis$parameters=="log(psi/(1-psi))",2:4])/(1+exp(boot_cis[boot_cis$parameters=="log(psi/(1-psi))",2:4]))
  boot_cis[boot_cis$parameters=="log(psi/(1-psi))",1] <- "psi"
  boot_cis[boot_cis$parameters=="log(omega2_nugg)",2:4] <- exp(boot_cis[boot_cis$parameters=="log(omega2_nugg)",2:4])
  boot_cis[boot_cis$parameters=="log(omega2_nugg)",1] <- "omega2_nugg"
  
  # xi
  boot_cis[str_detect(boot_cis$parameters,"xi"),2:4] <- exp(boot_cis[str_detect(boot_cis$parameters,"xi"),2:4])
  
  return(boot_cis)
}

get.scaling.value.fn <- function(control, rat, human, boot_cis){
  sp <- function(x,k) max(0,x-k)
  sp <- Vectorize(sp)
  
  # read in controls
  multi.xi.on <- control$multi.xi.on[1]
  xi.var <- control$xi.var[1]
  cov.rat <- control$rat[!is.na(control$rat)]
  cov.rat.sub <- cov.rat[str_detect(cov.rat, "sp.", negate=TRUE) & str_detect(cov.rat, "as.factor", negate=TRUE)]
  cov.rat.sub <- unique(c(cov.rat.sub, unlist(str_extract_all(cov.rat[str_detect(cov.rat, ",", negate=TRUE)],  "(?<=\\().+?(?=\\))"))))
  cov.human <- control$human[!is.na(control$human)]
  if(multi.xi.on == TRUE){cov.human <- c(cov.human, paste0("as.factor(",xi.var,")"))}
  cov.human.sub <- cov.human[str_detect(cov.human, "sp.", negate=TRUE) & str_detect(cov.human, "as.factor", negate=TRUE)]
  cov.human.sub <- unique(c(cov.human.sub, unlist(str_extract_all(cov.human[str_detect(cov.human, ",", negate=TRUE)],  "(?<=\\().+?(?=\\))")), cov.rat.sub))
  formula.rat <- as.formula(paste("~-1+",paste(cov.rat,collapse="+")))
  formula.human <- as.formula(paste("outcome ~",paste(cov.human,collapse="+")))
  
  rat <- na.omit(rat[,c("X","Y","data_type","outcome","offset","offset_req", cov.rat.sub)])
  human <- na.omit(human[,c("X","Y","outcome",cov.human.sub)])
  
  # rattiness covariate model matrix at all locations
  rat.human <- bind_rows(rat[,cov.rat.sub], human[,cov.rat.sub])
  D.aux <- as.matrix(model.matrix(formula.rat, data=rat.human))
  
  # coordinates
  coords.rat <- as.matrix(rat[,c("X","Y")])
  coords.human <- as.matrix(human[,c("X","Y")])
  coords.set <- data.frame(rbind(coords.rat,coords.human))
  colnames(coords.set)  <- c("X","Y")
  id.rat <- 1:nrow(coords.rat)
  id.human <- (nrow(coords.rat)+1):(nrow(coords.rat)+nrow(coords.human))
  ID <- create.ID.coords(coords.set,~X+Y)
  coords <- unique(coords.set)
  
  # standardise rattiness covariates
  p <- ncol(D.aux)
  N <- nrow(coords)
  D <- matrix(NA,nrow=N,ncol=p)
  for(i in 1:p) {
    D[,i] <- tapply(D.aux[,i],ID,max)
  }
  D.unscale <- D
  D <- scale(D)
  out1 <- attr(D,"scaled:scale")
  names(out1) <- colnames(D.aux)
  
  glm.fit <- glm(formula.human, data=human,family=binomial,x=TRUE)
  D.human.us <- glm.fit$x 
  p.human <- ncol(D.human.us)
  D.human.s <- scale(D.human.us[,2:p.human])
  out2 <- attr(D.human.s,"scaled:scale")
  return(c(out1,out2))
}
