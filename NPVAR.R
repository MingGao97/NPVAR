library("np")
library("mgcv")


#' @description Estimate topological ordering for nonparametric DAG with equal variance assumption
#'
#'
#' @param x A n x p data matrix/dataframe w/o names
#' @param run.np Use kernel regression to estimate conditional variance
#' @param run.mgcv Use GAM to estimate conditional variance
#'  If both are True, conduct both and output comparison for two methods
#'  Use GAM for final result
#' @param verbose Whether to output detailed estimation process
#' @param layer.select Whether to do layer selection
#' @param eta If layer.select is True, use this as fixed eta for layer selection
#' @param x2 Same format as but independent with x. If layer.select is True and eta is NULL,
#'  use x2 for estimation for eta adaptively
#' @param eta.scaler A scaler to control the scale of eta when variance of the system is large.
#'  The final eta = eta.scaler * eta_hat.
#'
#'
#' @return If layer.select is True, return a list contained with a vector ancestors of ordering,
#'  and a list of layers. If not, return a vector of ancestors of ordering.
NPVAR = function(x,
                 run.np = FALSE,
                 run.mgcv = TRUE,
                 verbose = FALSE,
                 layer.select = FALSE,
                 eta = NULL,
                 x2 = NULL,
                 eta.scaler = 1)
{
  if (layer.select & is.null(eta) & is.null(x2)) {
    stop("Dataset X2 shoud be provided if doing adaptive layer selection!")
  }

  if(!is.data.frame(x)){
    x = data.frame(x)
  }
  n = nrow(x)
  p = ncol(x)
  node.index = 1:p
  names(node.index) = names(x)

  ### Find the first source by minimum marginal variance
  condvars = apply(x, 2, var)
  source = argmin(condvars)
  if (layer.select) {
    layers = list()
    diff.condvars = condvars - min(condvars)
    if (is.null(eta)) {
      if(!is.data.frame(x2)){
        x2 = data.frame(x2)
      }
      condvars2 = apply(x2, 2, var)
      eta_t = max(abs(condvars - condvars2)) * eta.scaler
      layers[[1]] = node.index[diff.condvars <= eta_t]
    } else {
      layers[[1]] = node.index[diff.condvars <= eta]
    }
    ancestors = layers[[1]]
  } else {
    ancestors = c(source)
  }

  if (verbose) {
    print(condvars)
  }

  ### Main loop of algorithm
  l = 2
  while(length(ancestors) < p - 1){
    ### Estimate conditional variance for each descendant
    descendants = node.index[-ancestors]
    condvars = est.condvars(x, descendants, node.index, ancestors,
                            verbose, run.np, run.mgcv)

    ### Choose the node with minimum conditional variance as source
    min_index = argmin(condvars)
    source_next = node.index[names(condvars)[min_index]]
    if (layer.select) {
      diff.condvars = condvars - min(condvars)
      if (is.null(eta)) {
        condvars2 = est.condvars(x2, descendants, node.index, ancestors,
                                 verbose, run.np, run.mgcv)
        eta_t = max(abs(condvars - condvars2)) * eta.scaler
        layers[[l]] = node.index[names(condvars)[diff.condvars <= eta_t]]
      } else {
        layers[[l]] = node.index[names(condvars)[diff.condvars <= eta]]
      }
      ancestors = c(ancestors, layers[[l]])
    } else {
      ancestors = c(ancestors, source_next)
    }

    l = l + 1 # layer indicator

    if(verbose){
      print(condvars)
    }
  }

  ### Last inclusion
  descendants = (node.index)[-ancestors]
  if (length(ancestors) < p) {
    ancestors = c(ancestors, descendants)
    if (layer.select) {
      layers[[l]] = descendants
    }
  }

  ### Final output
  print(names(ancestors))
  if (layer.select) {
    return(list(ancestors = ancestors, layers = layers))
  } else {
    return(ancestors)
  }
}


### A helper function to estimate conditional variance for descendants
### Parameters are the same in main function
est.condvars = function(x, descendants, node.index, ancestors, verbose, run.np, run.mgcv) {

  condvars = rep(NA, length(descendants))
  names(condvars) = names(descendants)

  for(j in seq_along(descendants)){

    current_node = descendants[j]
    if (verbose){
      message(sprintf("Checking %s ~ %s", toupper(names(node.index)[current_node]),
                      paste(toupper(names(ancestors)), collapse =" + ")))
    }

    ### Compute nonparametric regression using kernel / np package
    if(run.np){
      npformula = "x[,current_node] ~ "
      for(a in ancestors){
        npformula = paste0(npformula, "x[,", a, "] + ")
      }
      npformula = substr(npformula, start = 1, stop = nchar(npformula) - 3)
      npformula = as.formula(npformula)

      bw.all = npregbw(formula = npformula,
                       regtype = "ll",
                       bwmethod = "cv.aic",
                       data = x)
      model.np = npreg(bws = bw.all)
      fit.np = predict(model.np,
                       data = x,
                       newdata = x
      )
      condvar.np = var(x[,current_node]) - var(fit.np)
    }

    ### Compute nonparametric regression using GAM / mgcv package
    if(run.mgcv){
      mgcvformula = "x[,current_node] ~ "
      for(a in ancestors){
        mgcvformula = paste0(mgcvformula, "s(x[,", a, "], bs='ps', sp=0.6) + ")
      }
      mgcvformula = substr(mgcvformula, start = 1, stop = nchar(mgcvformula) - 3)
      mgcvformula = as.formula(mgcvformula)

      b1 = mgcv::gam(mgcvformula, data = x)
      fit.gam = predict(b1)

      condvar.gam = var(x[,current_node]) - var(fit.gam)
    }

    if(run.mgcv && run.np){
      message(sprintf("np estimate: %f / gam estimate: %f / difference: %f",
                      condvar.np, condvar.gam, condvar.np - condvar.gam))
    }

    ### Record the estimated conditional variance
    if(run.np){
      condvars[j] = condvar.np
    } else if(run.mgcv){
      condvars[j] = condvar.gam
    } else{
      stop("Invalid nonparametric regression choice!")
    }
  }
  if(any(condvars < 0)) stop("Error: Negative variance found")

  return(condvars)
}

### Compute the argmin set of a vector
argmin = function(x) which(x == min(x))
