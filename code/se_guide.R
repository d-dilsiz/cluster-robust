# University of Basel
# Düzgün Dilsiz
# 2024-03-09

# -----------------------------------------------------------------------------
# Set working directory, load packages and data
# -----------------------------------------------------------------------------

# Working directory
# setwd("...")

# Packages
library(haven)
library(AER)
library(plm)

# Open data
rm(list=ls())
crime <- data(Crime)

# -----------------------------------------------------------------------------
# Model Crime
# -----------------------------------------------------------------------------

# fixed effects estimation with non-parametric time trend, plm
plm <- plm(polpc ~ crmrte + density + year,
           data=Crime, index=c("county","year"), model="within")

# define cluster vector with same length as model matrix (needed for later)
cluster <- Crime$county[complete.cases(Crime[ , c("density","crmrte")])]
complete <- complete.cases(Crime[,c("density","crmrte")])

# dimensions
NT <- length(plm$model[,1])
N <- length(unique(cluster))
T <- length(unique(Crime$year[complete]))
K <- length(plm$coefficients) # no. of time-varying variables (incl. T dummies)

# -----------------------------------------------------------------------------
# vcov = sss (cluster-robust) [STATA]
# -----------------------------------------------------------------------------

cluster_robust_se_sss_stata <- function(model, cluster_array){
  
  # finite data modification sss (fdm)
  fdm <-  (N / (N - 1)) * ((NT - 1) / (NT - K - 1))
  
  # take transformed residuals
  uhat <- residuals(model)
  
  # omega_i is Xi'üi'üiXi of cluster i
  omega_i <- lapply(unique(cluster_array), function(clusterid) {
    i <- cluster_array == clusterid
    
    X_i <- model.matrix(model)[i, , drop = FALSE] 
    # FALSE ensures dimensions are not dropped for clusters with one observation
    # time-invariant variables get a zero. For unbalanced panels, include 
    # complete to i.
    
    X_i <- X_i[ ,colSums(X_i != 0) > 0]
    # remove time-invariant variables: only then this command is needed, in case 
    # of unbalanced panel it deletes observations with just one observation
    
    t(X_i) %*% tcrossprod(uhat[i]) %*% X_i
  })
  
  # omega is sum over all clusters omega_i, times fdm
  omega <- fdm * Reduce('+', omega_i) # for meat
  
  # sandwich formula
  X <- model.matrix(model) # for bread
  X <- X[ ,colSums(X != 0) > 0] # remove time-invariant variables 
  vcov <- solve(t(X) %*% X, diag(K)) %*% omega %*% solve(t(X) %*% X, diag(K))
  
  # function returns vcov-matrix 
  return(vcov)
}
coeftest(plm, vcov=cluster_robust_se_sss_stata(plm, cluster))

# -----------------------------------------------------------------------------
# vcov = sss (cluster-robust) [R]
# -----------------------------------------------------------------------------

cluster_robust_se_sss <- function(model, cluster_array){
  
  # finite data modification sss (fdm)
  fdm <-  (N / (N - 1)) * ((NT - 1) / (NT - K))
  
  # take transformed residuals
  uhat <- residuals(model)
  
  # omega_i is Xi'üi'üiXi of cluster i
  omega_i <- lapply(unique(cluster_array), function(clusterid) {
    i <- cluster_array == clusterid
    
    X_i <- model.matrix(model)[i, , drop = FALSE] 
    
    X_i <- X_i[ ,colSums(X_i != 0) > 0]
    
    t(X_i) %*% tcrossprod(uhat[i]) %*% X_i
  })
  
  # omega is sum over all clusters omega_i, times fdm
  omega <- fdm * Reduce('+', omega_i)
   
  # sandwich formula
  X <- model.matrix(model)
  X <- X[ ,colSums(X != 0) > 0] 
  vcov <- solve(t(X) %*% X, diag(K)) %*% omega %*% solve(t(X) %*% X, diag(K))
  
  # function returns vcov-matrix 
  return(vcov)
}

coeftest(plm, vcov=cluster_robust_se_sss(plm, cluster))
coeftest(plm, vcov=vcovHC(plm, type="sss", cluster="group"))

# -----------------------------------------------------------------------------
# vcov = HC0 (cluster-robust)
# -----------------------------------------------------------------------------

cluster_robust_se_hc0 <- function(model, cluster_array){
  
  # finite data modification HC0 (fdm)
  fdm <-  1
  
  # take transformed residuals
  uhat <- residuals(model)
  
  # omega_i is Xi'üi'üiXi of cluster i
  omega_i <- lapply(unique(cluster_array), function(clusterid) {
    i <- cluster_array == clusterid
    
    X_i <- model.matrix(model)[i, , drop = FALSE] 
    
    X_i <- X_i[ ,colSums(X_i != 0) > 0]
    
    t(X_i) %*% tcrossprod(uhat[i]) %*% X_i
  })
  
  # omega is sum over all clusters omega_i, times fdm
  omega <- fdm * Reduce('+', omega_i) 
  
  # sandwich formula
  X <- model.matrix(model)
  X <- X[ ,colSums(X != 0) > 0] 
  vcov <- solve(t(X) %*% X, diag(K)) %*% omega %*% solve(t(X) %*% X, diag(K))
  
  # function returns vcov-matrix 
  return(vcov)
}

coeftest(plm, vcov=cluster_robust_se_hc0(plm, cluster))
coeftest(plm, vcov=vcovHC(plm, type="HC0", cluster="group"))

# -----------------------------------------------------------------------------
# vcov = HC1 (cluster-robust)
# -----------------------------------------------------------------------------

cluster_robust_se_hc1 <- function(model, cluster_array){
  
  # finite data modification HC0 (fdm)
  fdm <-  (NT / (NT - K)) # or equivalently: (N/(N-1)) * ((NT-T)/(NT-K))
  
  # take transformed residuals
  uhat <- residuals(model)
  
  # omega_i is Xi'üi'üiXi of cluster i
  omega_i <- lapply(unique(cluster_array), function(clusterid) {
    i <- cluster_array == clusterid
    
    X_i <- model.matrix(model)[i, , drop = FALSE] 
    
    X_i <- X_i[ ,colSums(X_i != 0) > 0]
    
    t(X_i) %*% tcrossprod(uhat[i]) %*% X_i
  })
  
  # omega is sum over all clusters omega_i, weighted by fdm
  omega <- fdm * Reduce('+', omega_i) 
  
  # sandwich formula
  X <- model.matrix(model) 
  X <- X[ ,colSums(X != 0) > 0] 
  vcov <- solve(t(X) %*% X, diag(K)) %*% omega %*% solve(t(X) %*% X, diag(K))
  
  # function returns vcov-matrix 
  return(vcov)
}

coeftest(plm, vcov=cluster_robust_se_hc1(plm, cluster))
coeftest(plm, vcov=vcovHC(plm, type="HC1", cluster="group"))

# -----------------------------------------------------------------------------
# vcov = HC2 (cluster-robust)
# -----------------------------------------------------------------------------

cluster_robust_se_hc2 <- function(model, cluster_array){
  
  # finite data modification HC2 (fdm)
  fdm <-  1
  
  # compute hat values
  X <- model.matrix(model)
  X <- X[ ,colSums(X != 0) > 0]
  hatvalues <- diag(X %*% solve(t(X) %*% X) %*% t(X))
  
  hat_mean <- mean(hatvalues)
  uhat <- resid(model) / sqrt((1 - hatvalues))
  
  # omega_i is Xi'üi'üiXi of cluster i
  omega_i <- lapply(unique(cluster_array), function(clusterid) {
    i <- cluster_array == clusterid
    
    X_i <- model.matrix(model)[i, , drop = FALSE] 
    
    X_i <- X_i[ ,colSums(X_i != 0) > 0]
    
    t(X_i) %*% tcrossprod(uhat[i]) %*% X_i
  })
  
  # omega is sum over all clusters omega_i, times fdm
  omega <- fdm * Reduce('+', omega_i)
  
  # sandwich formula
  X <- model.matrix(model) 
  X <- X[ ,colSums(X != 0) > 0] 
  vcov <- solve(t(X) %*% X, diag(K)) %*% omega %*% solve(t(X) %*% X, diag(K))
  
  # function returns vcov-matrix
  return(vcov)
}

coeftest(plm, vcov=cluster_robust_se_hc2(plm, cluster))
coeftest(plm, vcov=vcovHC(plm, type="HC2", cluster="group"))

# -----------------------------------------------------------------------------
# vcov = HC3 (cluster-robust)
# -----------------------------------------------------------------------------

cluster_robust_se_hc3 <- function(model, cluster_array){
  
  # finite data modification HC2 (fdm)
  fdm <-  1
  
  # compute hat values
  X <- model.matrix(model)
  X <- X[ ,colSums(X != 0) > 0]
  hatvalues <- diag(X %*% solve(t(X) %*% X) %*% t(X))

  uhat <- (resid(model)) / ((1 - hatvalues))
  
  # omega_i is Xi'üi'üiXi of cluster i
  omega_i <- lapply(unique(cluster_array), function(clusterid) {
    i <- cluster_array == clusterid
    
    X_i <- model.matrix(model)[i, , drop = FALSE] 
    
    X_i <- X_i[ ,colSums(X_i != 0) > 0]
    
    t(X_i) %*% tcrossprod(uhat[i]) %*% X_i
  })
  
  # omega is sum over all clusters omega_i, times fdm
  omega <- fdm * Reduce('+', omega_i) 
  
  # sandwich formula
  X <- model.matrix(model) 
  X <- X[ ,colSums(X != 0) > 0] 
  vcov <- solve(t(X) %*% X, diag(K)) %*% omega %*% solve(t(X) %*% X, diag(K))
  
  # function returns vcov-matrix 
  return(vcov)
}

coeftest(plm, vcov=cluster_robust_se_hc3(plm, cluster))
coeftest(plm, vcov=vcovHC(plm, type="HC3", cluster="group"))

