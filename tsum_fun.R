## Function to summarize data from multiple trials
##
## Using only base functionality for simplicity

## Input is a data.frame with a response and a grouping variable

tsum <- function(x, var.names = c("lrr","Trial_ID"), 
                 method = c("ici","lm","min-max"),
                 level = 0.95, order = TRUE,
                 add.id = TRUE){
  
  if(!class(x) == "data.frame") stop("only for data.frames")
  
  method <- match.arg(method)
  
  if(method == "ici"){
    ## Calcualte mean by trial
    frm <- as.formula(paste0(var.names[1],"~",var.names[2]))
    ans <- aggregate(formula = frm, data = x, FUN = mean)
    ans$lb <- aggregate(formula = frm, data = x, FUN = lb_fun, level = level)[,2]
    ans$ub <- aggregate(formula = frm, data = x, FUN = ub_fun, level = level)[,2]
    ans$n <- aggregate(formula = frm, data = x, FUN = length)[,2]
    names(ans) <- c(var.names[2], "m","lb","ub","n")
  }
  
  if(method == "lm"){
    ## Calcualte intervals using lm
    frm <- as.formula(paste0(var.names[1],"~",var.names[2]))
    nt <- aggregate(formula = frm, data = x, FUN = length)[,2]
    fit <- lm(formula = frm, data = x)
    ndat <- data.frame(unique(x[,var.names[2]]))
    names(ndat) <- var.names[2]
    prd <- predict(fit, newdata = ndat, interval = "confidence")
    ans <- cbind(ndat, prd)
    ans$n <- nt
    names(ans) <- c(var.names[2], "m","lb","ub","n")
  }
  
  if(method == "min-max"){
    ## Calcualte mean by trial
    rsp <- x[,var.names[1]]
    grp <- x[,var.names[2]]
    ans <- aggregate(rsp ~ grp, data = x, FUN = mean)
    ans$lb <- aggregate(rsp ~ grp, data = x, FUN = min)
    ans$ub <- aggregate(rsp ~ grp, data = x, FUN = max)
    ans$n <- aggregate(rsp ~ grp, data = x, FUN = length)
    names(ans) <- c(var.names[2], "m","lb","ub","n")
  }
  
  if(order & add.id){
    ans2 <- ans[order(ans[,"m"]),]
    ans2$id <- 1:nrow(ans2)
    return(ans2)
  }else{
    return(ans)
  }
}

lb_fun <- function(x, level = 0.95){
  n.x <- length(x)
  m.x <- mean(x)
  se.x <- sd(x)/sqrt(n.x)
  qnt <- 1 - (1 - level)/2
  tv <- qt(qnt, n.x-1)
  lbv <- m.x -  tv * se.x
  return(lbv)
}

ub_fun <- function(x, level = 0.95){
  n.x <- length(x)
  m.x <- mean(x)
  se.x <- sd(x)/sqrt(n.x)
  qnt <- 1 - (1 - level)/2
  tv <- qt(qnt, n.x-1)
  ubv <- m.x +  tv * se.x
  return(ubv)
}