## This function will calculate a prediction interval in the context of
## meta-analysis
## As defined in Higgins et al. (2009) it considers the parameteric uncertainty
## and the between study uncertainty

## Implemented cases:
##
## object = numeric
## object = lmer
## object = MCMCglmm
## object = data.frame
##
## Alternative 'boot' method for 'lmerMod'
## Methods available for 'MCMCglmm'
## - tdist
## - mcmc
## - simulate
## - predict

## This is version 0.2

pred_int <- function(x, level = 0.95, 
                     interval = c("prediction","confidence"),
                     method = c("tdist","boot","mcmc","simulate","predict"),
                     m.method = "REML",
                     nsim = 500,...){
  
  interval <- match.arg(interval)
  method <- match.arg(method)
  
  if(class(x) != "numeric" && 
     class(x) != "lmerMod" &&
     class(x) != "MCMCglmm" &&
     class(x) != "data.frame") stop("object not supported")
  
  if(class(x) == "lmerMod" && interval == "confidence"){ 
    stop("not implemented, use confint instead")
  }
  
  if(class(x) == "numeric"){
    ans <- predict(lm(x ~ 1), newdata=data.frame(x = mean(x)), 
                   interval = interval, level = level)
  }
  
  if(class(x) == "data.frame"){
    ans <- pred_int_df(x, method = m.method, 
                       var.names = NULL,
                       level = level,
                       interval = interval)
  }
  
  if(class(x) == "lmerMod" && method == "tdist"){
    ## Extract needed components
    n.k <- x@Gp[2] ## Number of trials
    alph <- (1 - level)/2
    t.val <- qt(1 - alph, n.k-2) ## t.value
    mu <- fixef(x)
    tau.sqrd <- VarCorr(x)[[1]][1] ## between studies variance
    se.mu <- summary(x)$coefficients[2] ## Intercept standard error
    ## The following line combines parametric uncertainty (se.mu)
    ## With sampling uncertainty (tau.sqrd)
    se.muh <- sqrt(tau.sqrd + se.mu^2) 
    pdi <- t.val * se.muh
    mu.lci <- mu - pdi
    mu.uci <- mu + pdi
    ans <- as.vector(c(mu, mu.lci, mu.uci))
  }
  
  if(class(x) == "lmerMod" && method == "boot"){
    tmp <- suppressWarnings(bootMer(x, 
                                    nsim = nsim, 
                                    pred_int, 
                                    use.u = TRUE)$t)
    ans <- colMeans(tmp)
  }
 
  if(class(x) == "MCMCglmm" && method == "tdist"){
    ans <- pred_int_mcg(x, level = level)
  }
  
  if(class(x) == "MCMCglmm" && method == "mcmc"){
    mu.chain <- as.vector(x$Sol[,1])
    tau.chain <- sqrt(as.vector(x$VCV[,1]))
    err.chain <- rnorm(length(mu.chain), 0, tau.chain)
    prd.chain <- mu.chain + err.chain
    mu.pi <- c(median(prd.chain), quantile(prd.chain, probs = c(0.025,0.975)))
    ans <- mu.pi
  }
  
  if(class(x) == "MCMCglmm" && method == "simulate"){
    ## Extract number of effective samples
    nrs <- nrow(as.matrix(x$Sol))
    ndat <- data.frame(lrr = rep(0,nrs), 
                       Trial_ID = "A_new_trial")
    prd.b <- scale(simulate(x, newdata = ndat,
                           marginal = NULL), scale = FALSE)
    pd.mu <- as.matrix(x$Sol[,1])
    prdi.c.mu <- pd.mu + prd.b
    ans <- quantile(prdi.c.mu, probs = c(0.5, 0.025,0.975))
  }
  
  if(class(x) == "MCMCglmm" && method == "predict"){
    ## I will assume the pr = TRUE when the model was ran
    ## It is assumed that there is a factor called 'Trial_ID'
    if(as.character(x$Random$formula)[2] != "Trial_ID")
      stop("Trial_ID not found")
    tns <- x$Z@Dimnames[[2]]
    trnms <- sapply(tns, FUN = function(x) strsplit(x, ".", fixed = TRUE)[[1]][2])
    rsp.nm <- as.character(x$Fixed$formula[2])
    ## The next line is the median for the intercept
    rsp.m <- as.data.frame(summary(soy.mc$Sol)[[2]][1,3])
    names(rsp.m) <- rsp.nm
    ndat <- data.frame(rsp.m, 
                       Trial_ID = as.vector(trnms))
    ans <- colMeans(predict(x, 
                            newdata = ndat, 
                            marginal = NULL, 
                            interval = interval))
  }
  
  names(ans) <- c("fit","lwr","upr")
  return(ans)
}

pred_int_mcg <- function(x, level = 0.95){
  ## This function extracts elements from a MCMCglmm object
  ## and calculates a prediction interval
  if(class(x) != "MCMCglmm") stop("only for MCMCglmm objects")
                                  
  n.k <- length(x$Z@Dimnames[[2]])
  alph <- (1 - level)/2
  t.val <- qt(1 - alph, n.k-2) ## t.value
  mu <- summary(x)$solution[1]
  pvm <- summary(x$VCV)[[1]]
  tau.sqrd <- pvm[1] ## between studies variance
  se.mu <- summary(x$Sol)[[1]][1,2] ## Intercept se
  se.muh <- sqrt(tau.sqrd + se.mu^2) 
  pdi <- t.val * se.muh
  mu.lci <- mu - pdi
  mu.uci <- mu + pdi
  ans <- as.vector(c(mu, mu.lci, mu.uci))
  return(ans)
}

pred_int_df <- function(x, method = "REML", 
                        interval = c("prediction", "confidence"),
                        var.names = c("TRT_Yld","CTR_Yld","Trial_ID"),
                        level = 0.95){
  
  interval <- match.arg(interval)
  ## This method assumes there is a data.frame
  ## with some names
  if(class(x) != "data.frame") stop("only for data frames")
  ## Make sure names align
  x2 <- x[,var.names]
  if(ncol(x2) != 3) stop("var.names might be wrong")
  ## Calculate preliminaries
  x.trt <- x[,var.names[1]]
  x.ctr <- x[,var.names[2]]
  trial <- x[,var.names[3]]
  x.m1i <- aggregate(x.trt ~ trial, data = x, FUN = mean)[,2]
  x.m2i <- aggregate(x.ctr ~ trial, data = x, FUN = mean)[,2]
  x.sd1i <- aggregate(x.trt ~ trial, data = x, FUN = sd)[,2]
  x.sd2i <- aggregate(x.ctr ~ trial, data = x, FUN = sd)[,2]
  x.n1i <- aggregate(x.trt ~ trial, data = x, FUN = length)[,2]
  x.n2i <- aggregate(x.ctr ~ trial, data = x, FUN = length)[,2]
  
  rr.escl <- escalc("ROM", m1i = x.m1i, m2i = x.m2i, 
                    sd1i = x.sd1i, sd2i = x.sd2i,
                    n1i = x.n1i, n2i = x.n2i)
  
  x.rma <- rma(yi, vi, data = rr.escl, method = method, weighted = FALSE)
  
  x.pred <- unclass(predict(x.rma, level = level))
  if(interval == "prediction"){
    ans <- c(x.pred$pred, x.pred$cr.lb, x.pred$cr.ub)
  }else{
    ans <- c(x.pred$pred, x.pred$ci.lb, x.pred$ci.ub)
  }
  names(ans) <- c("fit","lwr","upr")
  return(ans)
}
