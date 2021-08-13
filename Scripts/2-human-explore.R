# Script: 2-human-explore.R
# Authors: Max T Eyre, Emanuele Giorgi, Peter Diggle 
# Date: 05/05/2021
# Article: 
# Explanation: In this script we explore the human leptospiral infection data, providing summaries of key variables
# and exploring the functional form of their relationship with the log-odds of probability of infection. Univariable
# associations are described followed by a variable selection process to identify key variables of interest for the
# multivariable analysis. Then the interaction between rattiness and household elevation level on infection risk 
# is explored. Finally, we test for evidence of residual spatial correlation after controlling for variables.

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
library(interactions)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(MuMIn)

sp <- function(x,k) max(0,x-k)
sp <- Vectorize(sp)

## Step 1 - Load human data and summarise main characteristics ----
human <- read_csv("Data/human_data.csv") %>%
  mutate_at(c("sex","valley","elev_levels","hillside","open_sewer_10m",
              "exposed_to_sewer","work_constr", "work_salesman","work_trashcollect",
              "work_mud","work_floodwater","work_sewer","floodwater_freq_bin",
              "sew_water_freq_bin","mud_freq_bin"),as.factor)

## Step 2 - Explore functional form of continuous variables using GAMs ----
# age & income
gam.age.income <- gam(data=human, outcome~s(age) + as.factor(sex)   + 
                  as.factor(valley)  + s(hh_income_pcap_pess30_dailydol), 
                  family=binomial("logit"), method="REML")
gamviz.age.income <- getViz(gam.age.income) 
plot(gamviz.age.income, allTerms = T) # age: knot at 30yrs, income linear

# relative elevation
gam.elev <- gam(data=human, outcome~s(age) + as.factor(sex)   + 
                  as.factor(valley)  + s(hh_income_pcap_pess30_dailydol) + 
                  s(rel_elev), family=binomial("logit"), method="REML")
gamviz.elev <- getViz(gam.elev) 
plot(gamviz.elev, allTerms = T) # knot at 20m

# land cover
gam.lc <- gam(data=human, outcome~s(age) + as.factor(sex)   + 
                  as.factor(valley)  + s(hh_income_pcap_pess30_dailydol) + 
                s(lc_30m_imperv), family=binomial("logit"), method="REML")
gamviz.lc <- getViz(gam.lc)
plot(gamviz.lc, allTerms = T) # linear

# years of education
gam.educ <- gam(data=human, outcome~s(age) + as.factor(sex)   + 
                as.factor(valley)  + s(hh_income_pcap_pess30_dailydol) + 
                s(educ_yrs), family=binomial("logit"), method="REML")
gamviz.educ <- getViz(gam.educ)
plot(gamviz.educ, allTerms = T) # knot at 5yrs

## Step 3 - univariable analysis ----

OR <- sapply(c("sex","hh_income_pcap_pess30_dailydol","illiteracy","lc_30m_imperv","open_sewer_10m","exposed_to_sewer",
                "work_constr","work_salesman","work_trashcollect","work_mud","work_floodwater",
                "work_sewer","hillside"),
              
              function(var) {
                
                formula    <- as.formula(paste("outcome ~ ", var, "+ (1|ID)"))
                res.logist <- glmer(formula, data = human, family = binomial, nAGQ = 8, 
                                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
                n <- nrow(coef(summary(res.logist)))
                x <- confint(res.logist)
                
                out <- c(coef(summary(res.logist))[2:n,],x[(2+1):(n+1),1:2])
                return(out)
              })

# 3 level variables
var3 <- c("age + sp(age,30)","rel_elev + sp(rel_elev,20)","educ_yrs + sp(educ_yrs,5)","valley","elev_levels","floodwater_freq_bin","sew_water_freq_bin","mud_freq_bin")
var3out <- tibble()
for (i in 1:length(var3)){
  formula    <- as.formula(paste("outcome ~ ", var3[i], "+ (1|ID)"))
  res.logist <- glmer(formula, data = human, family = binomial, nAGQ = 8, 
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  x <- confint(res.logist)
  out <- cbind(coef(summary(res.logist))[2:3,],x[3:4,1:2])
  var3out <- rbind(var3out,out)
}

lc30 <- OR[,"lc_30m_imperv"]/10
exp(lc30)

OR.full <- rbind(t(OR),var3out) %>%
  mutate(Estimate = round(exp(Estimate),2), 
         `2.5 %` = round(exp(`2.5 %`),2), 
         `97.5 %` = round(exp(`97.5 %`),2)) %>%
  dplyr::select(Estimate, `2.5 %`, `97.5 %`)

# adjusted
aOR <- sapply(c("illiteracy","lc_30m_imperv","open_sewer_10m","exposed_to_sewer",
               "work_constr","work_salesman","work_trashcollect","work_mud","work_floodwater",
               "work_sewer","hillside"),
             
             function(var) {
               
               formula    <- as.formula(paste("outcome ~ age + sp(age,30) + sex + as.factor(valley)  + hh_income_pcap_pess30_dailydol +", var, "+ (1|ID)"))
               res.logist <- glmer(formula, data = human, family = binomial, nAGQ = 8, 
                                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
               n <- nrow(coef(summary(res.logist)))
               x <- confint(res.logist,method="Wald")
               
               out <- c(coef(summary(res.logist))[8,],x[9,1:2])
               return(out)
             })

# 3 level variables
var3aOR <- c("rel_elev + sp(rel_elev,20)","educ_yrs + sp(educ_yrs,5)","elev_levels","floodwater_freq_bin","sew_water_freq_bin","mud_freq_bin")
var3outaOR <- tibble()
for (i in 1:length(var3aOR)){
  formula    <- as.formula(paste("outcome ~ age + sp(age,30) + sex + as.factor(valley)  + hh_income_pcap_pess30_dailydol +", var3aOR[i], "+ (1|ID)"))
  res.logist <- glmer(formula, data = human, family = binomial, nAGQ = 8, 
                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
  x <- confint(res.logist,method="Wald")
  out <- cbind(coef(summary(res.logist))[8:9,],x[9:10,1:2])
  var3outaOR <- rbind(var3outaOR,out)
}

# get for confounders
res.logist <- glmer(outcome ~ age + sp(age,30) + sex + as.factor(valley)  + hh_income_pcap_pess30_dailydol + (1|ID), data = human, family = binomial, nAGQ = 8, 
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
x <- confint(res.logist,method="Wald")
out <- cbind(coef(summary(res.logist))[2:7,],x[3:8,1:2])

var3outaOR <- rbind(var3outaOR,out)

lc30aOR <- aOR[,"lc_30m_imperv"]/10
exp(lc30aOR)

aOR.full <- rbind(t(aOR),var3outaOR) %>%
  mutate(Estimate = round(exp(Estimate),2), 
         `2.5 %` = round(exp(`2.5 %`),2), 
         `97.5 %` = round(exp(`97.5 %`),2)) %>%
  dplyr::select(Estimate, `2.5 %`, `97.5 %`)


## Step 4 - Model selection ----
# a priori = age, sex, valley, income

### AIC SELECTION: within domains ----
# social status

human$age.sp30 <- sp(human$age, 30)
data.socialstatus <- human %>% filter(!is.na(educ_yrs))
model.socialstatus <- glmer(data =data.socialstatus, 
                            outcome~  age + sp(age,30) + as.factor(sex) + as.factor(valley) + hh_income_pcap_pess30_dailydol +
                              educ_yrs + sp(educ_yrs,5) + as.factor(illiteracy) + (1|ID), nAGQ = 8, 
                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), family=binomial)

options(na.action = "na.fail") 
d.socialstatus <- dredge(global.model = model.socialstatus, fixed = ~age + sp(age,30) + as.factor(sex) + as.factor(valley) + hh_income_pcap_pess30_dailydol)
options(na.action = "na.omit")
# SELECT: none

# household environment
library(parallel)
library(doParallel)
data.hhenv <- human %>% filter(!is.na(hillside))
data.hhenv$rel_elev.sp20 <- sp(data.hhenv$rel_elev, 20)
cl <- makeCluster(6)
registerDoParallel(cl)
clusterExport(cl,"data.hhenv")
clusterEvalQ(cl, library(lme4))

model.hhenv <- glmer(data=data.hhenv, 
                     outcome~  age + age.sp30 + as.factor(sex) + as.factor(valley) + hh_income_pcap_pess30_dailydol +
                       lc_30m_imperv + as.factor(elev_levels) +  hillside + rel_elev + rel_elev.sp20 +
                       open_sewer_10m + exposed_to_sewer + (1|ID), nAGQ = 8, 
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), family=binomial())
options(na.action = "na.fail") 
d.hhenv <- pdredge(global.model = model.hhenv, 
                   fixed = ~age + age.sp30 + as.factor(sex) + as.factor(valley) + hh_income_pcap_pess30_dailydol,
                   cluster = cl,
                   subset=  (rel_elev  & rel_elev.sp20) | !rel_elev.sp20)
options(na.action = "na.omit")
# SELECT: lc_30m_imperv

# occupational risk
data.occuprisk <- human
cl <- makeCluster(6)
registerDoParallel(cl)
clusterExport(cl,"data.occuprisk")
clusterEvalQ(cl, library(lme4))

model.occuprisk <- glmer(data=data.occuprisk, 
                     outcome~  age + age.sp30 + as.factor(sex) + as.factor(valley) + hh_income_pcap_pess30_dailydol +
                       work_constr + work_salesman+ work_trashcollect+ work_mud+ 
                       work_floodwater+ work_sewer + (1|ID), nAGQ = 8, 
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), family=binomial())
options(na.action = "na.fail") 
d.occuprisk <- pdredge(global.model = model.occuprisk, 
                       fixed = ~age + age.sp30 + as.factor(sex) + as.factor(valley) + hh_income_pcap_pess30_dailydol,
                       cluster = cl)
options(na.action = "na.omit")
# SELECT: work_salesman (delta == 0)

# behavioural exposure
data.expose <- human %>% filter(!is.na(floodwater_freq_bin),!is.na(sew_water_freq_bin))
cl <- makeCluster(5)
registerDoParallel(cl)
clusterExport(cl,"data.expose")
clusterEvalQ(cl, library(lme4))

model.expose <- glmer(data=data.expose, 
                         outcome~  age + age.sp30 + as.factor(sex) + as.factor(valley) + hh_income_pcap_pess30_dailydol +
                        as.factor(floodwater_freq_bin) + as.factor(sew_water_freq_bin) +  (1|ID), nAGQ = 8, 
                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), family=binomial())
options(na.action = "na.fail") 
d.expose <- pdredge(global.model = model.expose, fixed = ~age + age.sp30 + as.factor(sex) + as.factor(valley) + hh_income_pcap_pess30_dailydol,
                    cluster = cl)
options(na.action = "na.omit")
# SELECT: floodwater_freq_bin

# put together all selected variables and reduce
data.all <- human %>% filter(!is.na(floodwater_freq_bin),!is.na(work_salesman))
cl <- makeCluster(6)
registerDoParallel(cl)
clusterExport(cl,"data.all")
clusterEvalQ(cl, library(lme4))

model.all <- glmer(data=data.all, 
                   outcome ~ age + age.sp30 + as.factor(sex) + as.factor(valley) + hh_income_pcap_pess30_dailydol +
                  lc_30m_imperv + work_salesman + floodwater_freq_bin +  
                    elev_levels*pred.mean.rattiness + (1|ID),  family=binomial,
                   nAGQ = 8, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
options(na.action = "na.fail") 
d.all <- pdredge(global.model = model.all, fixed = ~age + age.sp30 + as.factor(sex) + as.factor(valley) + hh_income_pcap_pess30_dailydol +
                   elev_levels*pred.mean.rattiness,
                    cluster = cl)
# SELECT: work_salesman + floodwater_freq_bin

# final model
human$elev_levels <- relevel(human$elev_levels,ref = c("3"))
model.human.mixed <- glmer(data=human, 
                           outcome ~ age + sp(age, 30) + as.factor(sex) + as.factor(valley) + 
                             hh_income_pcap_pess30_dailydol + work_salesman +  floodwater_freq_bin +  
                             elev_levels*pred.mean.rattiness + (1|ID), family=binomial, nAGQ = 8, 
                           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(model.human.mixed)

## Step 5 - explore interaction with rattiness ----
human$elev_levels <- relevel(human$elev_levels,ref = c("1"))
levels(human$elev_levels)
model.human.interact <- glmer(data=human, 
                           outcome ~ age + sp(age, 30) + as.factor(sex) + as.factor(valley) + 
                             hh_income_pcap_pess30_dailydol + work_salesman +  floodwater_freq_bin +  
                             elev_levels*pred.mean.rattiness + (1|ID), family=binomial, nAGQ = 8, 
                           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))

p.inter <- interact_plot(model.human.interact, pred = pred.mean.rattiness, line.thickness = 0.8,
                         modx = elev_levels, interval = TRUE, int.width = 0.95, outcome.scale ="link", int.type = "confidence", 
                         legend.main = NULL, pred.labels = c("Rattiness"), 
                         modx.labels = c("Low (0 - 6.7m)", "Medium (6.7 - 15.6m)", "High (>15.6m)"), colors = c("#56B4E9","#56B4E9","#56B4E9"), vary.lty = FALSE) + 
  theme_bw() + xlab("Rattiness") + ylab("Log-odds infection risk") +
  theme(text = element_text(size=txtsize), legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) + 
  scale_x_continuous(breaks = c(-1,0,1), labels = c(-1,0,1)) + 
  facet_wrap(~elev_levels)


# Step 6 - test for residual spatial correlation ----
human$units.m <- 1
human <- as.data.frame(human)
spatcor <- spat.corr.diagnostic(data= human, 
                     outcome~  age + sp(age, 30) + as.factor(sex) + as.factor(valley) + 
                       hh_income_pcap_pess30_dailydol + as.factor(work_salesman) +  as.factor(floodwater_freq_bin) +  
                       as.factor(elev_levels),
                     coords=~X+Y,
                     ID.coords = create.ID.coords(human, coords= ~X+Y),
                     likelihood="Binomial",
                     uvec=seq(0,100,length=20),
                     n.sim = 10000,
                     nAGQ = 4)
# no evidence of spatial correlation

