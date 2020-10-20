source('ANM_gp.R')
library(igraph)

### Simulated data generation
### Given graph type, source variance, number of nodes, sample size, graph degree, x2
### Return a list consists of data matrix and true graph adjacency matrix
data_simu = function(graph_type, errvar, d, n, s0, x2 = F){
  if (graph_type == 'MC-SIN') {
    G = markov_chain(d)
    X = sampleFromSin(G, n, errvar)
    if (x2) X2 = sampleFromSin(G, n, errvar)
  } else if (graph_type == 'MC-GP') {
    G = markov_chain(d)
    if (x2) {
      X = sampleDataFromG(2 * n, G, errvar = errvar)
      X2 = X[(n+1):(2*n),]
      X = X[1:n,]
    } else {
      X = sampleDataFromG(n, G, errvar = errvar)
    }
  } else if (graph_type == 'ER-AGP') {
    if ((d==5) & (s0>1)) {
      G = as.matrix(sparsebnUtils::random.graph(d, 9))
    } else {
      G = as.matrix(sparsebnUtils::random.graph(d, s0*d))
    }
    if (x2) {
      X = sampleDataFromG(2 * n, G, errvar = errvar)
      X2 = X[(n+1):(2*n),]
      X = X[1:n,]
    } else {
      X = sampleDataFromG(n, G, errvar = errvar)
    }
  } else if (graph_type == 'ER-SIN') {
    if ((d==5) & (s0>1)) {
      G = as.matrix(sparsebnUtils::random.graph(d, 9))
    } else {
      G = as.matrix(sparsebnUtils::random.graph(d, s0*d))
    }
    X = sampleFromSin(G, n, errvar)
    if (x2) X2 = sampleFromSin(G, n, errvar)
  } else if (graph_type == 'ER-NGP') {
    if ((d==5) & (s0>1)) {
      G = as.matrix(sparsebnUtils::random.graph(d, 9))
    } else {
      G = as.matrix(sparsebnUtils::random.graph(d, s0*d))
    }
    if (x2) {
      X = sampleDataFromG(2*n, G, errvar = errvar, parsFuncType=list(B=randomB(G),kap=0.01,sigmax=1,sigmay=1,output=FALSE))    
      X2 = X[(n+1):(2*n),]
      X = X[1:n,]
    } else {
      X = sampleDataFromG(n, G, errvar = errvar, parsFuncType=list(B=randomB(G),kap=0.01,sigmax=1,sigmay=1,output=FALSE))    
    }
  } else if (graph_type == 'SF-AGP') {
    G = as_adjacency_matrix(sample_pa(d, m = s0),sparse = F)
    if (x2) {
      X = sampleDataFromG(2 * n, G, errvar = errvar)
      X2 = X[(n+1):(2*n),]
      X = X[1:n,]
    } else {
      X = sampleDataFromG(n, G, errvar = errvar)
    }
  } else if (graph_type == 'SF-SIN') {
    G = as_adjacency_matrix(sample_pa(d, m = s0),sparse = F)
    X = sampleFromSin(G, n, errvar)
    if (x2) X2 = sampleFromSin(G, n, errvar)
  } else if (graph_type == 'SF-NGP') {
    G = as_adjacency_matrix(sample_pa(d, m = s0),sparse = F)
    if (x2) {
      X = sampleDataFromG(2*n, G, errvar = errvar, parsFuncType=list(B=randomB(G),kap=0.01,sigmax=1,sigmay=1,output=FALSE))    
      X2 = X[(n+1):(2*n),]
      X = X[1:n,]
    } else {
      X = sampleDataFromG(n, G, errvar = errvar, parsFuncType=list(B=randomB(G),kap=0.01,sigmax=1,sigmay=1,output=FALSE))    
    }
  }
  if (x2) {
    return(list(X=X, G=G, X2=X2))
  } else {
    return(list(X=X, G=G))
  }
}

### Transform some output into CPDAG
toCPDAG = function(est){
  before.zero = which(est == 0)
  after.zero = which((est - t(est)) == 0)
  cp = setdiff(after.zero, before.zero)  
  est[cp] = -1
  return(est)
}

### Estimate adjacency matrix using estimated ordering and significance given by GAM
prune = function(x, est_order, cutoff = 0.001){
  if(!is.data.frame(x)){
    x = data.frame(x)
  }
  p = dim(x)[2]
  adj = matrix(0, p, p)
  for (i in 2:p) {
    node = est_order[i]
    ancestors = est_order[1:(i-1)]
    
    mgcvformula = "x[,node] ~ "
    for(a in ancestors){
      mgcvformula = paste0(mgcvformula, "s(x[,", a, "], bs='ps', sp=0.6) + ")
    }
    mgcvformula = substr(mgcvformula, start = 1, stop = nchar(mgcvformula) - 3)
    mgcvformula = as.formula(mgcvformula)
    
    mod = mgcv::gam(mgcvformula)
    parents = ancestors[summary(mod)$s.pv < cutoff]
    adj[parents, node] = 1
  }
  return(adj)
}

### Calculate SHD given adjacency matrix
SHD = function(G_est, G_true){
  pred = which(G_est == 1)
  cond = which(G_true != 0)
  cond_reversed = which(t(G_true) != 0)
  extra = setdiff(pred, cond)
  reverse = intersect(extra, cond_reversed)
  gig_true = G_true + t(G_true)
  gig_est = G_est + t(G_est)
  gig_true[upper.tri(gig_true)] = 0
  gig_est[upper.tri(gig_est)] = 0
  pred_lower = which(gig_true != 0)
  cond_lower = which(gig_est != 0)
  extra_lower = setdiff(pred_lower, cond_lower)
  missing_lower = setdiff(cond_lower, pred_lower)
  shd = length(extra_lower) + length(missing_lower) + length(reverse)
  return(shd)
}

### Generate random Markov chain given number of nodes
markov_chain = function(d){
  adj = diag(d-1)
  adj = cbind(0,adj)
  adj = rbind(adj,0)
  random_sort = sample(1:d)
  adj = adj[random_sort, random_sort]
  return(adj)
}

### Compute the ordering from adjacency matrix
computeCausOrder <- function(G)
  # Copyright (c) 2013  Jonas Peters  [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms.
{
  p <- dim(G)[2]
  remaining <- 1:p
  causOrder <- rep(NA,p)
  for(i in 1:(p-1))
  {
    root <- min(which(colSums(G) == 0))
    causOrder[i] <- remaining[root]
    remaining <- remaining[-root]
    G <- G[-root,-root]
  }
  causOrder[p] <- remaining[1]
  return(causOrder)
}

### Sample from SIN model given adjacency matrix, sample size, source variance
sampleFromSin = function(G, n, errvar = 0.5){
  p = dim(G)[2]
  X = matrix(NA,n,p)
  causOrder = computeCausOrder(G)
  for (node in causOrder) {
    paOfNode = which(G[,node] == 1)
    if(length(paOfNode) == 0){
      X[,node] = rnorm(n, 0, sqrt(errvar))
    }else if(length(paOfNode) == 1){
      X[,node] = sin(X[,paOfNode]) + rnorm(n, 0, sqrt(errvar))
    }else{
      X[,node] = apply(sin(X[,paOfNode]), 1, sum) + rnorm(n, 0, sqrt(errvar))
    }
  }
  return(X)
}

### Sample from f=x*sin(x) model given adjacency matrix, sample size, source variance
sampleFromXSin = function(G, n, errvar = 0.5){
  p = dim(G)[2]
  X = matrix(NA,n,p)
  causOrder = computeCausOrder(G)
  for (node in causOrder) {
    paOfNode = which(G[,node] == 1)
    if(length(paOfNode) == 0){
      X[,node] = rnorm(n, 0, sqrt(errvar))
    }else if(length(paOfNode) == 1){
      X[,node] = X[,paOfNode]*sin(X[,paOfNode]) + rnorm(n, 0, sqrt(errvar))
    }else{
      X[,node] = apply(X[,paOfNode]*sin(X[,paOfNode]), 1, sum) + rnorm(n, 0, sqrt(errvar))
    }
  }
  return(X)
}

### Sample from f=x^2 model given adjacency matrix, sample size, source variance
sampleFromXsquare = function(G, n, errvar = 0.5){
  p = dim(G)[2]
  X = matrix(NA,n,p)
  causOrder = computeCausOrder(G)
  for (node in causOrder) {
    paOfNode = which(G[,node] == 1)
    if(length(paOfNode) == 0){
      X[,node] = rnorm(n, 0, sqrt(errvar))
    }else if(length(paOfNode) == 1){
      X[,node] = X[,paOfNode]^2 + rnorm(n, 0, sqrt(errvar))
    }else{
      X[,node] = apply(X[,paOfNode]^2, 1, sum) + rnorm(n, 0, sqrt(errvar))
    }
  }
  return(X)
}

### Sample from f=x^{1/2} model given adjacency matrix, sample size, source variance
sampleFromSquareRoot = function(G, n, errvar = 0.5){
  p = dim(G)[2]
  X = matrix(NA,n,p)
  causOrder = computeCausOrder(G)
  for (node in causOrder) {
    paOfNode = which(G[,node] == 1)
    if(length(paOfNode) == 0){
      X[,node] = rnorm(n, 0, sqrt(errvar))
    }else if(length(paOfNode) == 1){
      X[,node] = sqrt(abs(X[,paOfNode])) + rnorm(n, 0, sqrt(errvar))
    }else{
      X[,node] = apply(sqrt(abs(X[,paOfNode])), 1, sum) + rnorm(n, 0, sqrt(errvar))
    }
  }
  return(X)
}

### Test if there is violation of estimated ordering under the true adjacency matrix
### Return whether ordering is correct and count of violations
test_order = function(topo_order, adj){
  edges = which(adj == 1, arr.ind = T)
  order_Index = order(topo_order)
  count = 0
  for (i in 1:nrow(edges)) {
    if (order_Index[edges[i,1]] > order_Index[edges[i,2]]) {
      count = count + 1
    }
  }
  return(list(right = as.numeric(count==0), count = count))
}

### Equal Variance algorithm
### Code from EqVarDAG from Wenyu Chen
EqVarDAG_TD_internal<-function(X){
  n<-dim(X)[1]
  p<-dim(X)[2]
  done<-NULL
  done<-p+1
  S<-cov(X)
  Sinv<-solve(S)
  for(i in 1:p){
    varmap<-seq(p)[-done]
    v<-which.min(diag(solve(Sinv[-done,-done])))[1]
    done<-c(done,varmap[v])
  }
  return(list(TO=done[-1],support=NULL))
}