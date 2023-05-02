rm(list = ls())
ipm_data=readRDS("C:/Users/Dave/OneDrive - University of Wyoming/Traits/WD_Traits/shapefile_work/gn_data.RDS")
#### Lines 5-38 are data manipulation unique to DHA DCL project
unique(ipm_data$Site)
ipm_data$size_t = ipm_data$basalArea_genet

ipm_data$Site = gsub("Ft Valley - COC-S1A", "FVEF S1", ipm_data$Site)
ipm_data$Site = gsub("Ft Valley - COC-S1B", "FVEF S1", ipm_data$Site)
ipm_data$Site = gsub("Ft Valley - COC-S2A", "FVEF S2", ipm_data$Site)
ipm_data$Site = gsub("Ft Valley - COC-S2B", "FVEF S2", ipm_data$Site)
ipm_data$Site = gsub("Ft Valley - COC-S3A", "FVEF S3", ipm_data$Site)

unique(ipm_data$Site)

ipm_data$code = ipm_data$species

ipm_data$code = gsub("Bouteloua gracilis", "BOUGRA", ipm_data$code)
ipm_data$code = gsub("Bromus ciliatus", "BROCIL", ipm_data$code)
ipm_data$code = gsub("Carex geophila", "CARGEO", ipm_data$code)
ipm_data$code = gsub("Elymus elymoides", "ELYELY", ipm_data$code)
ipm_data$code = gsub("Festuca arizonica", "FESARI", ipm_data$code)
ipm_data$code = gsub("Koeleria macrantha", "KOEMAC", ipm_data$code)
ipm_data$code = gsub("Muhlenbergia montana", "MUHMON", ipm_data$code)
ipm_data$code = gsub("Muhlenbergia tricholepis", "MUHTRI", ipm_data$code)
ipm_data$code = gsub("Poa fendleriana", "POAFEN", ipm_data$code)
ipm_data$code = gsub("Sporobolus interruptus", "SPOINT", ipm_data$code)

ipm_data$log_sizet = log(ipm_data$size_t)
ipm_data$log_sizetplus1 = log(ipm_data$size_tplus1)

saveRDS(ipm_data, "ipm_data.RDS")

##### GOOD FROM HERE
ipm_data = readRDS("ipm_data.RDS")

## new dataset to trim
ipm_data_trim = ipm_data
ipm_data_trim$survives_tplus1 = as.factor(ipm_data_trim$survives_tplus1)
## remove very large and small individuals to get rid of line patterns
#ipm_data_trim = ipm_data_trim[log(ipm_data_trim$size_t) > -11 & log(ipm_data_trim$size_tplus1) > -11,] ## REMOVES ALL mortality
#ipm_data_trim = ipm_data_trim[!is.na(ipm_data_trim$z_Year),]
## remove individuals growing or shrinking by a factor of CUTOFF
#cutoff = 5
#ipm_data_trim = ipm_data_trim[ipm_data_trim$size_t/ipm_data_trim$size_tplus1 < cutoff,] 
#ipm_data_trim = ipm_data_trim[ipm_data_trim$size_tplus1/ipm_data_trim$size_t < cutoff,]

#ipm_data_trim$prop_neighbor_cover_20cm = ipm_data_trim$neighbors_area_20cm/ipm_data_trim$nBuff_area_20cm
#ipm_data_trim$asin_prop_cover_20cm = asin(sqrt(ipm_data_trim$prop_neighbor_cover_20cm))
plot(log(ipm_data_trim$size_t), log(ipm_data_trim$size_tplus1), main = "Original")

## Removing extreme suspect points
test = ipm_data_trim[!(!log(ipm_data_trim$size_t) > -10 & !log(ipm_data_trim$size_tplus1)/log(ipm_data_trim$size_t) > .7),]
test = test[!(!log(test$size_tplus1) > -10 & !log(test$size_tplus1)/log(test$size_t) < 1.5),]
test = test[log(test$size_t) > -12.8 & log(test$size_tplus1) > -12.8,]

plot(log(test$size_t), log(test$size_tplus1), main = "Extreme suspects removed")
par(mfrow = c(1,1))
ipm_data_trim2 = test
ipm_data_trim2 = ipm_data_trim2[!grepl("NA", rownames(ipm_data_trim2)),]
ipm_split = split(ipm_data_trim2, list(ipm_data_trim2$species, ipm_data_trim2$Site, ipm_data_trim2$z_Year))


ipm_split = split(ipm_data_trim, list(ipm_data_trim$species, ipm_data_trim$Site, ipm_data_trim$z_Year))

ipm_split = ipm_split[lapply(ipm_split, nrow)>20] ## Filter down to only years with > 10 observations
#### READY TO GO

#############################################
### Fecundity estimates - flowering and seeds
##############################################
fec = read.csv(paste0(getwd(),"/Fecundity Data/fecundity_allSpp.csv"))
fec$log_area_m2 = log((fec$basal.area.cm2/10000))
fec = fec[!is.na(fec$log_area_m2),]
glume = read.csv(paste0(getwd(),"/Fecundity Data/Glumes.csv"))
glume$log_area_m2 = log((glume$BA.cm2/10000))

fec$flwr = ifelse(fec$numb.stalks>0, 1, 0)

fec_split = split(fec, fec$species)
glume_split = split(glume, glume$spp.code)

flower.reg = glm(flwr ~ log_area_m2, data=fec, family = "binomial")
flwr_fun = function(x){
  flower.reg = glm(flwr ~ log_area_m2, data=x, family = "binomial")
  out = list(flwr.int = coefficients(flower.reg)[1], flwr.slope = coefficients(flower.reg)[2])
  return(out)
}

flwr.params = lapply(fec_split, flwr_fun)

seeds_fun = function(x){
  seed.mod = glm(sum.glumes ~ log_area_m2, data=x, family = "poisson")
  return(seed.mod)
}

seed.mods = lapply(glume_split, seeds_fun)

########################################
#### Density independent, min N = 10
########################################
ipm_split = split(ipm_data_trim, list(ipm_data_trim$species, ipm_data_trim$Site, ipm_data_trim$z_Year))
ipm_split = ipm_split[lapply(ipm_split, nrow)>10] ## Filter down to only years with > 10 observations
paramList_n10 = list()
i = 0
j=i+1
for(i in j:length(ipm_split)){
  dat_i = ipm_split[[i]]
  #if(unique(dat_i$z_Year == "2021")){
  #  next
  #}
  sp.code = unique(dat_i$code)
  naames = stringr::str_split(names(ipm_split)[[i]], "[.]")
  seed.mod = seed.mods[[paste(sp.code)]]
  newdat = data.frame(log_area_m2 = dat_i$log_sizet)
  dat_i$seed.pred = round(predict(seed.mod, newdata = newdat, type = "response"),0)
  params=data.frame(
    #####
    species=naames[[1]][1],
    code = sp.code,
    site = naames[[1]][2],
    year = naames[[1]][3],
    #####
    surv.int=NA,
    surv.slope=NA,
    #####
    growth.int=NA,
    growth.slope=NA,
    growth.sd=NA,
    #####
    flwr.int=flwr.params[sp.code][[sp.code]]$flwr.int,
    flwr.slope=flwr.params[sp.code][[sp.code]]$flwr.slope,
    #####
    seed.int=coef(seed.mod)[1],
    seed.slope=coef(seed.mod)[2],
    #####
    recruit.size.mean=NA,
    recruit.size.sd=NA,
    #####
    establishment.prob=NA
  )
  
  surv.reg=glm(survives_tplus1 ~ log_sizet, data=dat_i,family= binomial, singular.ok = T)
  params$surv.int=coefficients(surv.reg)[1]
  params$surv.slope=coefficients(surv.reg)[2]
  
  growth.reg=lm(log_sizetplus1 ~ log_sizet, data=dat_i)
  params$growth.int=coefficients(growth.reg)[1]
  params$growth.slope=coefficients(growth.reg)[2]
  params$growth.sd=sd(resid(growth.reg))
  
  params$recruit.size.mean=mean(dat_i$log_sizet[dat_i$recruit==1], na.rm =TRUE)
  params$recruit.size.sd=sd(dat_i$log_sizet[dat_i$recruit==1], na.rm =TRUE)
  
  params$establishment.prob=sum(dat_i$recruit, na.rm = TRUE)/sum(dat_i$seed.pred,na.rm=TRUE)
  paramList_n10[[i]] = params

}

paramList_n10 = paramList_n10[lengths(paramList_n10)!=0]

########################################
#### Density independent, min N = 20
########################################
ipm_split = split(ipm_data_trim, list(ipm_data_trim$species, ipm_data_trim$Site, ipm_data_trim$z_Year))
ipm_split = ipm_split[lapply(ipm_split, nrow)>20] ## Filter down to only years with > 20 observations
paramList_n20 = list()
i = 0
j=i+1
for(i in j:length(ipm_split)){
  dat_i = ipm_split[[i]]
  #if(unique(dat_i$z_Year == "2021")){
  #  next
  #}
  sp.code = unique(dat_i$code)
  naames = stringr::str_split(names(ipm_split)[[i]], "[.]")
  seed.mod = seed.mods[[paste(sp.code)]]
  newdat = data.frame(log_area_m2 = dat_i$log_sizet)
  dat_i$seed.pred = round(predict(seed.mod, newdata = newdat, type = "response"),0)
  params=data.frame(
    #####
    species=naames[[1]][1],
    code = sp.code,
    site = naames[[1]][2],
    year = naames[[1]][3],
    #####
    surv.int=NA,
    surv.slope=NA,
    #####
    growth.int=NA,
    growth.slope=NA,
    growth.sd=NA,
    #####
    flwr.int=flwr.params[sp.code][[sp.code]]$flwr.int,
    flwr.slope=flwr.params[sp.code][[sp.code]]$flwr.slope,
    #####
    seed.int=coef(seed.mod)[1],
    seed.slope=coef(seed.mod)[2],
    #####
    recruit.size.mean=NA,
    recruit.size.sd=NA,
    #####
    establishment.prob=NA
  )
  
  surv.reg=glm(survives_tplus1 ~ log_sizet, data=dat_i,family= binomial, singular.ok = T)
  params$surv.int=coefficients(surv.reg)[1]
  params$surv.slope=coefficients(surv.reg)[2]
  
  growth.reg=lm(log_sizetplus1 ~ log_sizet, data=dat_i)
  params$growth.int=coefficients(growth.reg)[1]
  params$growth.slope=coefficients(growth.reg)[2]
  params$growth.sd=sd(resid(growth.reg))
  
  params$recruit.size.mean=mean(dat_i$log_sizet[dat_i$recruit==1], na.rm =TRUE)
  params$recruit.size.sd=sd(dat_i$log_sizet[dat_i$recruit==1], na.rm =TRUE)
  
  params$establishment.prob=sum(dat_i$recruit, na.rm = TRUE)/sum(dat_i$seed.pred,na.rm=TRUE)
  paramList_n20[[i]] = params
  
}

paramList_n20 = paramList_n20[lengths(paramList_n20)!=0]

##################################
## Density Dependent min N = 10
##################################
ipm_split = split(ipm_data_trim, list(ipm_data_trim$species, ipm_data_trim$Site, ipm_data_trim$z_Year))
ipm_split = ipm_split[lapply(ipm_split, nrow)>10] ## Filter down to only years with > 10 observations
paramList_dendep_n10 = list()

i = 0
j = i+1
for(i in j:length(ipm_split)){
  dat_i = ipm_split[[i]]
  if(unique(dat_i$z_Year == "2021")){
    next
  }
  sp.code = unique(dat_i$code)
  naames = stringr::str_split(names(ipm_split)[[i]], "[.]")
  seed.mod = seed.mods[[paste(sp.code)]]
  newdat = data.frame(log_area_m2 = dat_i$log_sizet)
  dat_i$seed.pred = round(predict(seed.mod, newdata = newdat, type = "response"),0)
  
  params_dd=data.frame(
    #####
    species=naames[[1]][1],
    code = sp.code,
    site = naames[[1]][2],
    year = naames[[1]][3],
    #####
    surv.int=NA,
    surv.slope=NA,
    surv.dd_beta = NA,
    surv.dd_p = NA,
    #####
    growth.int=NA,
    growth.slope=NA,
    growth.sd=NA,
    growth.dd_beta=NA,
    growth.dd_p = NA,
    #####
    flwr.int=flwr.params[sp.code][[sp.code]]$flwr.int,
    flwr.slope=flwr.params[sp.code][[sp.code]]$flwr.slope,
    #####
    seed.int=coef(seed.mod)[1],
    seed.slope=coef(seed.mod)[2],
    #####
    recruit.size.mean=NA,
    recruit.size.sd=NA,
    #####
    establishment.prob=NA
  )
  
  surv.reg=glm(survives_tplus1 ~ log_sizet + log(self_neighbors_area_20cm+0.00000000000001), data=dat_i, family = binomial, singular.ok = T)
  params_dd$surv.int=coefficients(surv.reg)[1]
  params_dd$surv.slope=coefficients(surv.reg)[2]
  params_dd$surv.dd_beta = coefficients(surv.reg)[3]
  params_dd$surv.dd_p = summary(surv.reg)$coefficients[3,4]
  
  growth.reg=lm(log_sizetplus1 ~ log_sizet + log(self_neighbors_area_20cm+0.00000000000001), data=dat_i)
  params_dd$growth.int=coefficients(growth.reg)[1]
  params_dd$growth.slope=coefficients(growth.reg)[2]
  params_dd$growth.sd=sd(resid(growth.reg))
  params_dd$growth.dd_beta = coefficients(growth.reg)[3]
  params_dd$growth.dd_p = summary(growth.reg)$coefficients[3,4]
  
  params_dd$recruit.size.mean=mean(dat_i$log_sizet[dat_i$recruit==1], na.rm =TRUE)
  params_dd$recruit.size.sd=sd(dat_i$log_sizet[dat_i$recruit==1], na.rm =TRUE)
  
  params_dd$establishment.prob=sum(dat_i$recruit, na.rm = TRUE)/sum(dat_i$seed.pred,na.rm=TRUE)
  paramList_dendep_n10[[i]] = params_dd
  
}
paramList_dendep_n10 = paramList_dendep_n10[lengths(paramList_dendep_n10)!=0]

##################################
## Density Dependent min N = 20
##################################
ipm_split = split(ipm_data_trim, list(ipm_data_trim$species, ipm_data_trim$Site, ipm_data_trim$z_Year))
ipm_split = ipm_split[lapply(ipm_split, nrow)>20] ## Filter down to only years with > 20 observations
paramList_dendep_n20_2 = list()

i = 0
j = i+1
for(i in j:length(ipm_split)){
  dat_i = ipm_split[[i]]
  if(unique(dat_i$z_Year == "2021")){
    next
  }
  sp.code = unique(dat_i$code)
  naames = stringr::str_split(names(ipm_split)[[i]], "[.]")
  seed.mod = seed.mods[[paste(sp.code)]]
  newdat = data.frame(log_area_m2 = dat_i$log_sizet)
  dat_i$seed.pred = round(predict(seed.mod, newdata = newdat, type = "response"),0)
  
  params_dd=data.frame(
    #####
    species=naames[[1]][1],
    code = sp.code,
    site = naames[[1]][2],
    year = naames[[1]][3],
    #####
    surv.int=NA,
    surv.slope=NA,
    surv.dd_beta = NA,
    surv.dd_p = NA,
    #####
    growth.int=NA,
    growth.slope=NA,
    growth.sd=NA,
    growth.dd_beta=NA,
    growth.dd_p = NA,
    #####
    flwr.int=flwr.params[sp.code][[sp.code]]$flwr.int,
    flwr.slope=flwr.params[sp.code][[sp.code]]$flwr.slope,
    #####
    seed.int=coef(seed.mod)[1],
    seed.slope=coef(seed.mod)[2],
    #####
    recruit.size.mean=NA,
    recruit.size.sd=NA,
    #####
    establishment.prob=NA
  )
  
  surv.reg=glm(survives_tplus1 ~ log_sizet + log(self_neighbors_area_20cm+0.00000000000001), data=dat_i, family = binomial, singular.ok = T)
  params_dd$surv.int=coefficients(surv.reg)[1]
  params_dd$surv.slope=coefficients(surv.reg)[2]
  params_dd$surv.dd_beta = coefficients(surv.reg)[3]
  params_dd$surv.dd_p = summary(surv.reg)$coefficients[3,4]
  
  growth.reg=lm(log_sizetplus1 ~ log_sizet + log(self_neighbors_area_20cm+0.00000000000001), data=dat_i)
  params_dd$growth.int=coefficients(growth.reg)[1]
  params_dd$growth.slope=coefficients(growth.reg)[2]
  params_dd$growth.sd=sd(resid(growth.reg))
  params_dd$growth.dd_beta = coefficients(growth.reg)[3]
  params_dd$growth.dd_p = summary(growth.reg)$coefficients[3,4]
  
  params_dd$recruit.size.mean=mean(dat_i$log_sizet[dat_i$recruit==1], na.rm =TRUE)
  params_dd$recruit.size.sd=sd(dat_i$log_sizet[dat_i$recruit==1], na.rm =TRUE)
  
  params_dd$establishment.prob=sum(dat_i$recruit, na.rm = TRUE)/sum(dat_i$seed.pred,na.rm=TRUE)
  paramList_dendep_n20_2[[i]] = params_dd
  
}
paramList_dendep_n20 = paramList_dendep_n20[lengths(paramList_dendep_n20)!=0]



L= 0.9*min(ipm_data_trim["log_sizet"], na.rm = T)
U= 1.1*max(ipm_data_trim["log_sizetplus1"], na.rm = T)
## number of cells in the discretized kernel
m=100
## boundary points (the edges of the cells defining the kernel)
b=L+c(0:m)*(U-L)/m 
## mesh points (midpoints of the cells)
y=0.5*(b[1:m]+b[2:(m+1)])
## width of the cells
h=diff(y[1:2])

###########################################
## Calculating lambda values from param lists
###########################################

###########################################
## Density Independent, min N  = 10
###########################################

i=0
j = i+1
for(i in j:length(paramList_n10)){
  params = paramList_n10[[i]]
  P <- h * (outer(y, y, P_z1z, params = params))
  F <- h * (outer(y, y, F_z1z, params = params))
  K <- P + F
  K[is.na(K)] = 0
  ## compute the eigen vectors / values 
  IPM.eig.sys <- eigen(K)
  ## lambda
  lambda <- Re(IPM.eig.sys$values[1])
  paramList_n10[[i]]["lambda"] = lambda
  
}
i==length(paramList_n10)

paramList_n10 = paramList_n10[lapply(paramList_n10, length) == 17]
temp = do.call(rbind.data.frame, paramList_n10)
lambdas_n10 = temp[complete.cases(temp),]

###########################################
## Density Independent, min N  = 20
###########################################

i=0
j = i+1
for(i in j:length(paramList_n20)){
  params = paramList_n20[[i]]
  P <- h * (outer(y, y, P_z1z, params = params))
  F <- h * (outer(y, y, F_z1z, params = params))
  K <- P + F
  K[is.na(K)] = 0
  ## compute the eigen vectors / values 
  IPM.eig.sys <- eigen(K)
  ## lambda
  lambda <- Re(IPM.eig.sys$values[1])
  paramList_n20[[i]]["lambda"] = lambda
  
}
i==length(paramList_n20)

paramList_n20 = paramList_n20[lapply(paramList_n20, length) == 17]
temp = do.call(rbind.data.frame, paramList_n20)
lambdas_n20 = temp[complete.cases(temp),]

###########################################
## Density Dependent, min N  = 10
###########################################

i=0
j = i+1
for(i in j:length(paramList_dendep_n10)){
  params = paramList_dendep_n10[[i]]
  P <- h * (outer(y, y, P_z1z, params = params))
  F <- h * (outer(y, y, F_z1z, params = params))
  K <- P + F
  K[is.na(K)] = 0
  ## compute the eigen vectors / values 
  IPM.eig.sys <- eigen(K)
  ## lambda
  lambda <- Re(IPM.eig.sys$values[1])
  paramList_dendep_n10[[i]]["lambda"] = lambda
  
}
i==length(paramList_dendep_n10)

paramList_dendep_n10 = paramList_dendep_n10[lapply(paramList_dendep_n10, length) == 21]
temp = do.call(rbind.data.frame, paramList_dendep_n10)
lambdas_dd_n10 = temp[complete.cases(temp),]

###########################################
## Density Dependent, min N  = 20
###########################################
i = 0
j = i+1
for(i in j:length(ipm_split)){
  dat_i = ipm_split[[i]]
  if(unique(dat_i$z_Year == "2021")){
    next
  }
  sp.code = unique(dat_i$code)
  naames = stringr::str_split(names(ipm_split)[[i]], "[.]")
  seed.mod = seed.mods[[paste(sp.code)]]
  newdat = data.frame(log_area_m2 = dat_i$log_sizet)
  dat_i$seed.pred = round(predict(seed.mod, newdata = newdat, type = "response"),0)
  
  params_dd=data.frame(
    #####
    species=naames[[1]][1],
    code = sp.code,
    site = naames[[1]][2],
    year = naames[[1]][3],
    #####
    surv.int=NA,
    surv.slope=NA,
    surv.dd_beta = NA,
    surv.dd_p = NA,
    #####
    growth.int=NA,
    growth.slope=NA,
    growth.sd=NA,
    growth.dd_beta=NA,
    growth.dd_p = NA,
    #####
    flwr.int=flwr.params[sp.code][[sp.code]]$flwr.int,
    flwr.slope=flwr.params[sp.code][[sp.code]]$flwr.slope,
    #####
    seed.int=coef(seed.mod)[1],
    seed.slope=coef(seed.mod)[2],
    #####
    recruit.size.mean=NA,
    recruit.size.sd=NA,
    #####
    establishment.prob=NA
  )
  
  surv.reg=glm(survives_tplus1 ~ log_sizet + log(self_neighbors_area_20cm+0.00000000000001), data=dat_i, family = binomial, singular.ok = T)
  params_dd$surv.int=coefficients(surv.reg)[1]
  params_dd$surv.slope=coefficients(surv.reg)[2]
  params_dd$surv.dd_beta = coefficients(surv.reg)[3]
  params_dd$surv.dd_p = summary(surv.reg)$coefficients[3,4]
  
  growth.reg=lm(log_sizetplus1 ~ log_sizet + log(self_neighbors_area_20cm+0.00000000000001), data=dat_i)
  params_dd$growth.int=coefficients(growth.reg)[1]
  params_dd$growth.slope=coefficients(growth.reg)[2]
  params_dd$growth.sd=sd(resid(growth.reg))
  params_dd$growth.dd_beta = coefficients(growth.reg)[3]
  params_dd$growth.dd_p = summary(growth.reg)$coefficients[3,4]
  
  params_dd$recruit.size.mean=mean(dat_i$log_sizet[dat_i$recruit==1], na.rm =TRUE)
  params_dd$recruit.size.sd=sd(dat_i$log_sizet[dat_i$recruit==1], na.rm =TRUE)
  
  params_dd$establishment.prob=sum(dat_i$recruit, na.rm = TRUE)/sum(dat_i$seed.pred,na.rm=TRUE)
  paramList_dendep_n20_2[[i]] = params_dd
  
}
paramList_dendep_n20 = paramList_dendep_n20[lengths(paramList_dendep_n20)!=0]

i=0
j = i+1
for(i in j:length(paramList_dendep_n20_2)){
  params = paramList_dendep_n20_2[[i]]
  P <- h * (outer(y, y, P_z1z, params = params))
  F <- h * (outer(y, y, F_z1z, params = params))
  K <- P + F
  K[is.na(K)] = 0
  ## compute the eigen vectors / values 
  IPM.eig.sys <- eigen(K)
  ## lambda
  lambda <- Re(IPM.eig.sys$values[1])
  paramList_dendep_n20_2[[i]]["lambda"] = lambda
  
}
i == length(paramList_dendep_n20)
paramList_dendep_n20_2 = paramList_dendep_n20_2[lapply(paramList_dendep_n20_2, length) == 21]
temp = do.call(rbind.data.frame, paramList_dendep_n20_2)
lambdas_dd_n20_2 = temp[complete.cases(temp),]

################################
## Merging into a single DF
################################

lams = data.frame(species = lambdas_n10$species,code = lambdas_n10$code, site = lambdas_n10$site, year = lambdas_n10$year, lam_n10 = lambdas_n10$lambda)
test = merge(x=lams, y=lambdas_n20[, c("species", "site", "year", "lambda")], by = c("species", "site", "year"), all.x = TRUE)
colnames(test)[6] = "lam_n20"
lams = test
test = merge(x=lams, y=lambdas_dd_n10[, c("species", "site", "year", "lambda")], by = c("species", "site", "year"), all.x = TRUE)
colnames(test)[7] = "lam_dd_n10"
lams = test
test = merge(x=lams, y=lambdas_dd_n20[, c("species", "site", "year", "lambda")], by = c("species", "site", "year"), all.x = TRUE)
colnames(test)[8] = "lam_dd_n20"
lams = test
View(lams)

saveRDS(lams, "lams.RDS")

summary(lm(lam_n10 ~ lam_dd_n10, data = lams))

plot(lam_dd_n10 ~ lam_n10, data = lams)
abline(0,1)

summary(lm(lam_n20 ~ lam_dd_n20, data = lams))

plot(lam_dd_n20 ~ lam_n20, data = lams)
abline(0,1)

table(lams$lam_n10 == 0)
table(lams$lam_n20 == 0)
table(lams$lam_dd_n10 == 0)
table(lams$lam_dd_n20 == 0)

table(na.omit(lams$lam_n10 > lams$lam_dd_n10))
table(na.omit(lams$lam_n20 > lams$lam_dd_n20))
