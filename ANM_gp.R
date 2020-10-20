# Taken from codeANM by Jonas Peters
# Code to generate data from a GP-nonlinear additive model SEM
#  Nonlinearities are sampled from a GP
#  sampleDataFromG is the main function, see documentation below
#

library(MASS)
sampleDataFromG <- function(n,G,funcType="GAM", 
                            parsFuncType=list(B=randomB(G),kap=1,sigmax=1,sigmay=1,output=FALSE), 
                            noiseType="normalRandomVariances", 
                            errvar = 0.5,
                            parsNoise=list(noiseExp=1,varMin=1,varMax=2,noiseExpVarMin=2,noiseExpVarMax=4,bound=rep(1,dim(G)[2])))
# Copyright (c) 2013-2013  Jonas Peters  [peters@stat.math.ethz.ch]
#                          Jan Ernest    [ernest@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.
#
#
# Generates n samples according to structural equation models based on the DAG G
# with specified function class and given noise distribution. Uses Gaussian processes to sample the
# nonlinear functions.
#
# INPUT:
#   m           number of samples that should be simulated
#   G           adjacency matrix of the full DAG
#   funcType    parameter to choose between different function types. Default is "GAM" which simulates from
#               an additive model (i.e. each node is an additive function of its parents). "GP" samples from
#               a fully non-additive function, wheras "GAMGP" interpolates between GAM and GP with parameter kap.
#               If funcType="linear" then the function samples from a linear SEM.
#   parsFuncType
#       kap     interpolates between a general additive model and a fully nonparametric model.
#               kap == 1 --> GAM with additive noise, kap == 0 --> fully nonparametric with additive noise
#       sigmax
#       sigmay  parameters used to simulate the Gaussian kernels
#   noiseType   specifies the type of additive noise in the model. Default is "normalRandomVariances" which simulates
#               Gaussian noise with random variances.
#   parsNoise   list of parameters to modify the noise distribution.
#
# OUTPUT:
#   X           n x p matrix containing the n samples of the specified model.

{
  if(funcType == "linear")
  {
    sampleDataFromGLinear(n=n,G=G,parsFuncType,noiseType,parsNoise)
    
  } else if(funcType == "GAM")
  {
    parsFuncType$kap = 1
    sampleDataFromGAMGP(n=n,G=G,errvar,parsFuncType,noiseType,parsNoise)
    
  } else if(funcType == "GP")
  {
    parsFuncType$kap = 0
    sampleDataFromGAMGP(n=n,G=G,errvar,parsFuncType,noiseType,parsNoise)
    
  } else if(funcType == "GAMGP")
  {
    sampleDataFromGAMGP(n=n,G=G,errvar,parsFuncType,noiseType,parsNoise)
  } else if(funcType == "Sigmoid")
  {
    sampleDataFromMonotoneSigmoid(n=n,G=G,parsFuncType,noiseType,parsNoise)
  } else
  {
    stop('This function type does not exist!')
  }
}

sampleDataFromGAMGP <- function(n, G, errvar, parsFuncType, noiseType, parsNoise)
  # INPUTS:   n:  number of samples
  #           G:  adjacency matrix of Graph to simulate from
  #           kap linearly interpolates between GAM and fully nonparametric model. -- kap == 1 --> GAM with additive noise, kap == 0 --> fully nonparametric with additive noise
  #           noiseExp: exponent to model deviations from Gaussian noise -- noiseExp == 1 --> Gaussian noise
  # OUTPUTS:  X:      sampled data
  
{
  p <- dim(G)[2]
  X <- matrix(NA,n,p)
  # determine the causal Order which is needed for sampling
  causOrder <- computeCausOrder(G)
  
  if(parsFuncType$output)
  {
    show(causOrder)
  }
  
  # sample noise variances
  # noiseVar <- runif(p,parsNoise$varMin,parsNoise$varMax)
  
  # change default to equal variances
  # errvar <- 0.5
  noiseVar <- rep(errvar, p)
  
  # loop through each node according to the causal order
  for(node in causOrder)
  {
    if(parsFuncType$output)
    {
      cat("generating GP for node ", node, "\r")
    }
    paOfNode <- which(G[,node] == 1)
    # simulation of noise at source nodes
    if(length(paOfNode) ==0)
    {
      if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
      {
        ran <- rnorm(n)
        noisetmp <- (sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)
      } else
      {
        error("This noiseType is not implemented yet.")
      }
      X[,node] <- noisetmp
    } else
    {
      nuPa <- length(paOfNode)
      X[,node] <- rep(0,n)
      
      # If kap>0 there is an additive model component
      if(parsFuncType$kap>0)
      {
        for(pa in paOfNode)
        {
          
          # sample parameters for Gaussian process
          kernPa <- computeGaussKernel(X[,pa],parsFuncType$sigmay,parsFuncType$sigmax)
          fpa <- mvrnorm(1,rep(0,n),kernPa)
          X[,node] <- X[,node] + parsFuncType$kap * fpa
        }
      }
      
      # if kap<1 there is a non-additive model component
      if(parsFuncType$kap<1)
      {
        kernAllPa <- computeGaussKernel(X[,paOfNode],parsFuncType$sigmay,parsFuncType$sigmax)
        fAllPa <- mvrnorm(1,rep(0,n),kernAllPa)
        if(parsFuncType$output & (parsFuncType$kap==0))
        {
          ### INCLUDE ADEQUATE PLOTTING FUNCTION (MOREDIMENSIONAL PLOTS) ###
        }
        X[,node] <- X[,node] + (1-parsFuncType$kap)*fAllPa
      }
      
      # Additive noise
      if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
      {
        ran <- rnorm(n)
        noisetmp <- (sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)
      } else
      {
        error("This noiseType is not implemented yet.")
      }
      X[,node] <- X[,node] + noisetmp
    }
  }
  
  return(X)
}

randomB <- function(G,lB = 0.1,uB = 0.9,twoIntervals = 1)
  # if twoIntervals == TRUE, lB and uB should be positive
  # Copyright (c) 2012-2012  Jonas Peters [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms.
{
  numCoeff <- sum(G)
  B <- t(G)
  if(numCoeff ==1)
  {
    coeffs <- sample(c(-1,1),size=numCoeff,0.5)^(twoIntervals) * runif(1,lB,uB)
  }
  else
  {
    coeffs <- diag(sample(c(-1,1),size=numCoeff,0.5)^(twoIntervals)) %*% runif(numCoeff,lB,uB)
  }
  B[B==1] <- coeffs
  return(B)
}

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

computeGaussKernel <- function(x, sigmay, sigmax)
{
  if(is.matrix(x)==FALSE){
    x<-as.matrix(x)}
  n <- dim(x)[1]
  
  xnorm<-as.matrix(dist(x,method="euclidean",diag=TRUE,upper=TRUE))
  xnorm<-xnorm^2
  
  KX <- sigmay * exp(-xnorm/(2*sigmax^2))
  
  return(KX)
}

