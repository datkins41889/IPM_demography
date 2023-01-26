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