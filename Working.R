ipm_data = GN_genet_05cm
ipm_data = ipm_data[,-c(6,7)]

ipm_data$neighbors_area_05cm = GN_genet_05cm$neighbors_area
ipm_data$nBuff_area_05cm = GN_genet_05cm$nBuff_area

ipm_data$neighbors_area_10cm = GN_genet_10cm$neighbors_area
ipm_data$nBuff_area_10cm = GN_genet_10cm$nBuff_area

ipm_data$neighbors_area_15cm = GN_genet_15cm$neighbors_area
ipm_data$nBuff_area_15cm = GN_genet_15cm$nBuff_area

ipm_data$neighbors_area_20cm = GN_genet_20cm$neighbors_area
ipm_data$nBuff_area_20cm = GN_genet_20cm$nBuff_area


## Build the discretized kernel
m=100
L= 0.9*min(dat$area_t)
U=1.1*max(dat$area_t)
mk_K <- function(m, params, L, U) {
  # mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  P <- h * (outer(meshpts, meshpts, P_z1z, params = params))
  F <- h * (outer(meshpts, meshpts, F_z1z, params = params))
  K <- P + F
  ## compute the eigen vectors / values 
  IPM.eig.sys <- eigen(K)
  ## lambda
  lambda <- Re(IPM.eig.sys$values[1])
  return(list(lambda=lambda, K = K, meshpts = meshpts, P = P, F = F))
}
mk_K(m=100,params = params, L = 0.9*min(dat$area_t), U=1.1*max(dat$area_t))