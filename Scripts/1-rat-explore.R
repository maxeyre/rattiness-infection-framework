# Script: 1-rat-explore.R
# Authors: Max T Eyre, Emanuele Giorgi, Peter Diggle 
# Date: 05/05/2021
# Article: 
# Explanation: In this script we explore the rat data by fitting the non-spatial rattiness model without covariates, 
# predicting Uj and exploring its relationship with covariates. Covariates are selected and then entered into the
# full spatial rattiness model. Finally predictions for rattiness are made at all household locations.

set.seed(180121)

library(geoR)
library(PrevMap)
library(randtoolbox)
library(numDeriv)
library(cowplot)
library(lme4)
library(mgcv)
library(tidyverse)
library(mgcViz)
library(paletteer)
library(pROC)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

### RAT DATA ----
source("Functions/1-rattiness_nonspatial_model_fns.R")

# Data in ----
rat <- read.csv("Data/rat_data.csv")

ID <- create.ID.coords(rat,~X+Y)
U.halton <- qnorm(halton(1000,1))

inst <- c("traps","plates","burrows","faeces","trails")
## ordered as 1 = traps, 2 = plates, 3 = burrows, 4 = faeces, 5 = trails
n.inst <- 5 # 5 indices
n.par <- 2*n.inst
par.start <- rep(0,n.par) # first guess at parameters

## Step 1 - Fit model (no covariates) and predict Uj ----
estim <- nlminb(par.start, 
                function(x) -llik.non.spatial.rattiness.fn(x),
                control=list(trace=1))
par0 <- estim$par
estim <- nlminb(par0, 
                function(x) -llik.non.spatial.rattiness.fn(x),
                control=list(trace=1))

## Step 2 - Predict mean rattiness at each location ----
n.loc <- length(unique(ID))
U.pred.mean <- sapply(1:n.loc,function(j) compute.pred.mean.Uj.no.cov.fn(estim$par,j))

# Join to covariate values
U.data <- data.frame(X=unique(rat[,c("X","Y")])[,1],
                     Y=unique(rat[,c("X","Y")])[,2],
                     U=U.pred.mean)

ratUj <- rat %>%
  dplyr::select(X, Y, valley, elevation, rel_elev, dist_sewer, dist_pri_sewer,
                dist_trash, lc_30m_veg, lc_30m_soil, lc_30m_imperv) %>%
  unique() %>%
  left_join(U.data, by=c("X","Y"))

## Step 3 - Explore covariate relationships with predicted rattiness using GAMs ----
# You can now plot Uj against covariates of interest (using ratUj) to explore relationships and decide how they 
# should be included in the model.
gam.rat <- gam(data=ratUj, U ~ s(rel_elev) + s(dist_trash)  + s(lc_30m_imperv))
gamviz.rat <- getViz(gam.rat)

width_ <- 120
height_ <- 100
res <- 300
txtsize <- 12

# relative elevation
jpeg("Outputs/Rat_GAM_rel_elev.jpeg", units="mm", width=width_, height=height_, res=res)
pp1 <- plot(gamviz.rat, select=1) + 
  l_ciPoly(alpha=0.7) +
  l_fitLine(linetype = 1)  +
  l_ciBar() +
  l_rug() +
  theme_bw() +
  xlab("Metres") +
  ylab("Rattiness") +
  theme(text = element_text(size=txtsize))
pp1
dev.off()

# distance to refuse
jpeg("Outputs/Rat_GAM_dist_trash.jpeg", units="mm", width=width_, height=height_, res=res)
pp2 <- plot(gamviz.rat, select=2) + 
  l_ciPoly(alpha=0.7) +
  l_fitLine(linetype = 1)  +
  l_ciBar() +
  l_rug() +
  theme_bw() +
  xlab("Metres") +
  ylab("Rattiness") +
  theme(text = element_text(size=txtsize))
pp2
dev.off()

# land cover
jpeg("Outputs/Rat_GAM_land_cover.jpeg", units="mm", width=width_, height=height_, res=res)
pp3 <- plot(gamviz.rat, select=3) + 
  l_ciPoly(alpha=0.7, fill="gray") +
  l_fitLine(linetype = 1)  +
  l_ciBar() +
  l_rug() +
  theme_bw() +
  xlab("Proportion impervious") +
  ylab("Rattiness") +
  theme(text = element_text(size=txtsize))
pp3
dev.off()

## Step 4 - Model selection using linear model ----
sp <- function(x,k) max(0,x-k)
sp <- Vectorize(sp)

model.rat <-  lm(data=ratUj, U~ rel_elev + sp(rel_elev,8) + sp(rel_elev,22) + dist_trash + sp(dist_trash, 50)  + lc_30m_imperv)

options(na.action = "na.fail") 
d.rat <- dredge(global.model = model.rat.final) %>%
  filter(!(is.na(`dist_trash`) & !is.na(`sp(dist_trash, 50)`)),
         !(is.na(`rel_elev`) & (!is.na(`sp(rel_elev, 8)`) | !is.na(`sp(rel_elev, 22)`))))
options(na.action = "na.omit")

## Step 5 - Fit non-spatial rattiness model with covariates to test for spatial correlation ----
D.aux <- as.matrix(model.matrix(~-1 + as.factor(valley) + rel_elev + sp(rel_elev,8) + sp(rel_elev,22) + dist_trash + sp(dist_trash, 50) + lc_30m_imperv, data=rat))

# standardise covariate values
coords <- unique(rat[,c("X","Y")])
p <- ncol(D.aux)
N <- nrow(coords)
D <- matrix(NA,nrow=N,ncol=p)
for(i in 1:p) {
  D[,i] <- tapply(D.aux[,i],ID,max)
}
D <- scale(D)

n.par <- 2*n.inst+p
par.start <- rep(0,n.par)

estim.full <- nlminb(par.start, 
                     function(x) -llik.non.spatial.rattiness.with.covariates.fn(x),
                     control=list(trace=1))

# Get empirical variogram when accounting for covariates
U.pred.mean <- sapply(1:n.loc,function(j) compute.pred.mean.Uj.with.cov.fn(estim.full$par,j))

U.data <- data.frame(X=unique(rat[,c("X","Y")])[,1],
                     Y=unique(rat[,c("X","Y")])[,2],
                     U=U.pred.mean)

ratUj_with_cov <- rat %>%
  dplyr::select(X, Y) %>%
  unique() %>%
  left_join(U.data, by=c("X","Y")) %>%
  as.data.frame(.)

vari.no.cov <- variogram(ratUj_with_cov,~U,~X+Y,uvec=seq(0,100,length=10))
env <- variog.mc.env(geodata=as.geodata(ratUj_with_cov),obj.variog=vari.no.cov,nsim=1000)

# Variogram
jpeg("Outputs/Rat_variogram.jpeg", units="mm", width=width_, height=height_, res=res)
pvar <- ggplot() + geom_ribbon(aes(x = env$u, ymin=env$v.lower, ymax=env$v.upper), alpha=0.7, fill="gray") +
  #geom_point(data = data.frame(u=vari.no.cov$u, v=vari.no.cov$v), aes(x=u, y=v)) +
  geom_line(data = data.frame(u=vari.no.cov$u, v=vari.no.cov$v), aes(x=u, y=v)) +
  theme_bw() +
  xlab("Distance (m)") +
  ylab("Semivariance") +
  theme(text = element_text(size=txtsize))
pvar
dev.off()

## Step 6 - Fit full rattiness spatial model ----
source(file = "Functions/2-rattiness_spatial_model_fn.R")
source(file = "Functions/scaling_fns.R")
model.name <- "R-rat_explore_R_H_spat_withnugg"
# Control sheet
control <- read_csv(paste0("Data/control_",model.name,".csv"))

output <- rattiness.eco.model.alt(control, rat,
                                  euclid.norm.method=TRUE, tol=0.5, n.iter.fit=10, rel.tol=1e-6,
                                  n.sim=c(50000), burnin=c(5000), 
                                  thin=10, MCMC.swap=FALSE, n.MCMC.swap=1)

## Step 7 - Predict rattiness at household locations----
# covariate values at household locations
pred.grid.human <- read_csv("Data/human_data.csv") %>%
  dplyr::select(ID, house_id, valley, X, Y, rel_elev, dist_trash, lc_30m_imperv)

source("Functions/3-rattiness_spatial_model_predict_fn.R")

crs.val <- CRS("+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs")

pred_out <- rattiness.ECO.predict(control, 
                                  rat, par_hat,
                                  n.sim = 55000, burnin = 5000, thin = 10,
                                  grid.pred, crs.val,
                                  get.raster = FALSE, n.pred.samples = 1)
pred.mean.rattiness <- pred_out$predictions.val$R.pred.hat
