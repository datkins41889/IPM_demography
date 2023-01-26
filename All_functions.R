## Helper functions for IPM script

calcParams <- function(
  dat = dat
){
  params=data.frame(
    surv.int=NA,
    surv.slope=NA,
    #####
    growth.int=NA,
    growth.slope=NA,
    growth.sd=NA,
    #####
    flwr.int=NA,
    flwr.slope=NA,
    #####
    seed.int=NA,
    seed.slope=NA,
    #####
    recruit.size.mean=NA,
    recruit.size.sd=NA,
    #####
    establishment.prob=NA
  )
  
  # 1. survival regression
  surv.reg=glm(survives_tplus1~area_t,data=dat,family=binomial)
  params$surv.int=coefficients(surv.reg)[1]
  params$surv.slope=coefficients(surv.reg)[2]
  
  # 2. growth regression
  growth.reg=lm(area_tplus1 ~ area_t,data=dat)
  params$growth.int=coefficients(growth.reg)[1]
  params$growth.slope=coefficients(growth.reg)[2]
  params$growth.sd=sd(resid(growth.reg))
  
  ## 3. flowering probability
  flower.reg = glm(flwr.sim ~ area_t, data=dat)
  params$flwr.int = coefficients(flower.reg)[1]
  params$flwr.slope = coefficients(flower.reg)[2]
  
  # 4. seeds regression
  seed.reg=glm(seed.sim ~ area_t, data = dat, family = "poisson")
  params$seed.int=coefficients(seed.reg)[1]
  params$seed.slope=coefficients(seed.reg)[2]
  
  
  # 5. size distribution of recruits
  # in the dataframe, recruits are those individuals who have a value for sizeNext but not for size
  params$recruit.size.mean=mean(dat$area_t[dat$recruit==1], na.rm =TRUE)
  params$recruit.size.sd=sd(dat$area_t[dat$recruit==1], na.rm =TRUE)
  
  ## 6. establishment probability
  ## 
  params$establishment.prob=sum(dat$recruit, na.rm = TRUE)/sum(dat$seed.sim,na.rm=TRUE)
  return(params)
}

calcElas = function(
  mk_Kout
) {
  IPM.sys = mk_Kout
  ## extract the mesh points and the kernel
  meshpts <- IPM.sys$meshpts
  h <- diff(meshpts[1:2])
  K <- IPM.sys$K
  P <- IPM.sys$P
  F <- IPM.sys$F
  
  
  ## compute the eigen vectors / values 
  IPM.eig.sys <- eigen(K)
  ## unnormalised stable size distribution ('w')
  w.z <- Re(IPM.eig.sys$vectors[,1])
  ## ... as before (nothing new here)
  lambda <- Re(IPM.eig.sys$values[1])
  
  ## unnormalised reproductive value distribution ('v')...
  v.z1 <- Re(eigen(t(K))$vectors[,1])
  
  ## kernel sensitivity
  K.sens <- outer(v.z1, w.z, "*")/sum(v.z1 * w.z * h)
  ## kernel elasticity (the divide-by-dz is needed to ensure we are working
  ## with the kernel, which is a function, and not the iteration matrix)
  K.elas <- K.sens * (K/h) / lambda
  ## does the kernel elasticity integrate to 1?...
  sum(K.elas) * h^2
  ## ...yes
  
  ## calculate the elasticities of the survival and reproduction components
  P <- IPM.sys$P / h
  F <- IPM.sys$F / h
  P.elas <- P * K.sens / lambda
  F.elas <- F * K.sens / lambda
  
  gs.elas = sum(P.elas) * h^2 ## ## Growth/Survival
  repro.elas = sum(F.elas) * h^2 ## Reproduction
  out = list(gs.elas, repro.elas)
  names(out) = c('Elas - GS', 'Elas - R')
  return(out)
  
}

############################################
g_z1z <- function(area_tplus1, area_t, params = params)
{
  mu <- params$growth.int + params$growth.slope * area_t           # mean size next year
  sig <- params$growth.sd                                    # sd about mean
  p.den.grow <- dnorm(area_tplus1, mean = mu, sd = sig)             # pdf that you are size area_tplus1 given you were size z
  return(p.den.grow)
}

## Survival function, logistic regression

s_z <- function(area_t, params)
{
  linear.p <- params$surv.int + params$surv.slope * area_t # linear predictor
  p <- 1/(1+exp(-linear.p))                            # logistic transformation to probability
  return(p)
}

## Probability of flowering function, logistic regression

p_bz <- function(area_t, params)
{
  linear.p <- params$flwr.int + params$flwr.slope * area_t      # linear predictor
  p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
  return(p)
}

## Recruitment function (N.B - from birth in spring to first summer), logistic regression
pr_z <- function(params) {
  linear.p <- params$recruit.size.mean                             # linear predictor
  p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
  return(p)
}

## Seed production function

b_z <- function(area_t, params)
{
  N <- exp(params$seed.int + params$seed.slope * area_t)    # seed production of a size z plant
  return(N)
}

## Recruit size pdf

c_0z1 <- function(area_t, params)
{
  mu <- params$recruit.size.mean
  sig <- params$recruit.size.sd
  p.deRecr <- dnorm(area_t, mean = mu, sd = sig)              # pdf of a size z1 recruit
  return(p.deRecr)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
P_z1z <- function (area_tplus1, area_t, params) {
  return( s_z(area_t, params) * g_z1z(area_tplus1, area_t, params) )
}

## Define the reproduction kernel
F_z1z <- function (area_tplus1, area_t, params) {
  
  return( p_bz(area_t, params) * b_z(area_t, params) * params$establishment.prob * c_0z1(area_tplus1, params))
  
}

## Quasi-Poisson draw
rqpois = function(n, mu, theta){
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}
simFlowers <- function(
  data=dat, fecund = fecund
){
  
  n <- nrow(data)
  data["flwr.sim"] = rbinom(n, 1, data["flwr.pred"][,1])
  data["glum.sim"] <- rqpois(n, data["glum.pred"][,1], sigma(mod))
  # simulate seeds per flower
  data["seed.sim"] = data["glum.sim"] * seedsPerGlume * data["flwr.sim"]
  return(data)
}



