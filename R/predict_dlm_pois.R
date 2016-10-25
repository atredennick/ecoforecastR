##' @name predict_dlm_pois
##' @title predict_dlm_pois
##' @export 
##' @author Mike Dietze and Andrew Tredennick
##' @description make predictions from fit dlm for new data
##' @param fit     fit_dlm output list
##' @param newdata array of driver data organized with dimensions [ensemble,time,variable]
##' @param n.iter  number of samples
##' @param steps   number of steps
##' @param start.time  time for start of forecast. Defaults to end of fit
##' @param include which source of uncertainty to include (vector of strings). Options are
##' \itemize{
##'  \item{I}{Initial Conditions}
##'  \item{P}{Parameters}
##'  \item{D}{Drivers}
##'  \item{E}{Process Error}
##' }
##' @return list
##' \itemize{
##'  \item{predict}{matrix of predictions [ens,time]}
##'  \item{index}{data.frame of the parameter (P) and driver (D) indices used for each ensemble member}
##' }
predict_dlm_pois <- function(fit,newdata=NULL,n.iter=5000,steps=NULL,start.time=NULL,include=c("I","P","D","E")){
  ## checks
  if(n.iter < 1){
    print("n.iter must be > 1")
    return(NULL)
  }

  
  ## set up variables
  if(is.null(steps)){
    if(is.null(newdata)){
      print("either newdata or steps needs to be provided")
      return(NULL)
    } else {
      steps = dim(newdata)[2]
    }
  }
  if(!("D" %in% include)){
    my.dims = dim(newdata)
    my.dims[1] = 1
    my.dimnames <- dimnames(newdata)
    my.dimnames[[1]] = 1
    newdata <- array(apply(newdata,3,apply,2,mean),
                        dim=my.dims,dimnames = my.dimnames)
  }
  
  params = as.matrix(fit$params)
  if(!("P" %in% include)){
    params <- as.matrix(apply(params,2,median))
  }
  if(ncol(params)==1) params <- t(params)
  
  IC = as.matrix(fit$predict)
  if(is.null(start.time)){
    start.time = ncol(IC)
  }
  IC = IC[,start.time]
  if(!("I" %in% include)){
    IC = median(IC)
  }
  
  ## set up storage
  predict = matrix(NA,n.iter,steps)
  
  ## sample indices from newdata
  index = data.frame(P = sample.int(nrow(params),n.iter,replace=TRUE),
                     D = sample.int(dim(newdata)[1],n.iter,replace=TRUE),
                     I = sample.int(length(IC),n.iter,replace=TRUE))
  x = IC[index$I]
  beta_IC = params[index$P,"beta_IC"]
  # beta    = params[index$P,paste0("beta",dimnames(newdata)[[3]])]
  
  if("E" %in% include){
    shape = params[index$P,"shape"]  ## convert from precision to SD
  } else {
    shape = params[index$P,"shape"]  ## convert from precision to SD
  }
  
  ## simulate
  for(t in 1:steps){
    Z  = newdata[index$D,t,]
    mu = beta_IC*x #+ apply( Z * beta,1,sum)
    x  = rgamma(n.iter, shape = shape, rate = shape/exp(mu))
    predict[,t] = rpois(n.iter,x)
  }

  ## output
  return(list(predict=predict,index=index,newdata=newdata))
}