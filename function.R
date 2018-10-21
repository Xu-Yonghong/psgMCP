#=============================================================================================================
# program: function.R
# note: code for psgMCP method, which has been proposed in manuscript "Identifying gene-environment interactions
#       incorporating prior information.
#
#=============================================================================================================


###############################################################################################################
#------------------------------------------------------------------------------------------------------------//
# Function: taylor.exp
# Input:                                                                                                     //
#           X                a matrix of dimension nobs*q, where q is the number of environment variables    //
#           W                a matrix of dimension nobs*(p*(q+1)), where p is the number of genes            //
#           Y                a numeric vector of lentgh n                                                    //
#           alpha0           a scaler respesenting the intercept                                             //
#           a0               a vector of length q, denoting the coefficient vector of X                      //
#           b0               a vector of length p*(q+1), representing the coefficient vector of W            //                                           //
#           type             response type.                                                                  //
# Output:                                                                                                    //
#    A list of length two: weight and Ystar                                                                  //
#           weight          weights for observations after Taylor expansion                                  //
#           Ystar           the new resonse after Taylor expansion                                           //
#------------------------------------------------------------------------------------------------------------//


#------------------------------------------------------------------------------------------------------------//
# Function: center.data
# Input:                                                                                                     //
#           X                a matrix of dimension nobs*q, where q is the number of environment variables    //
#           W                a matrix of dimension nobs*(p*(q+1)), where p is the number of genes            //
#           Y                a numeric vector of lentgh n                                                    //
#           alpha0           a scaler respesenting the intercept                                             //
#           a0               a vector of length q, denoting the coefficient vector of X                      //
#           b0               a vector of length p*(q+1), representing the coefficient vector of W            //                                          //
#           type             response type.                                                                  //
# Output:                                                                                                    //
#    A list of length four: Yols, Xols, Wols, meanlist                                                       //
#       - Yols               the centered response                                                           //
#       - Xols               the centered observations for E variables                                       //
#       - Wols               the centered observations for G and E-G variables                               //
#       - meanlist           a list of lenght three, including the weighted means of X, W,and Y                                                      //
#-------------------------------------------------------------------------------------------------------- ---//

#------------------------------------------------------------------------------------------------------------//
# Function: generate.data
# Input:                                                                                                     //
#           alpha0           a scaler respesenting the intercept                                             //
#           a                a vector of length q, denoting the coefficient vector of X                      //
#           b                a vector of length p*(q+1), representing the coefficient vector of W            //                                          //
#           n                the number of observations                                                      //
#           q                the number of E variables                                                       //
#           p                the number of G variables                                                       //
#           sigmaX           the correlation matrix for E variables, of dimension q*q                        //
#           sigmaZ           the correlation matrix for G variables, of dimension p*p                        //
#           sd.error         the standard deviation of error                                                 //
#           type             response type                                                                   //

# Output:                                                                                                    //
#           X                a matrix of dimension n*q                                                       //
#           W                a matrix of dimension n*(p*(q+1))                                               //
#           Z                a matrix of dimension n*p                                                       //
#           Y                a response vector of lentgh n                                                   //
#-------------------------------------------------=----------------------------------------------------------//

#------------------------------------------------------------------------------------------------------------//
# Function: cv.sgMCP.prior
# Input:                                                                                                     //
#           X                a matrix of dimension n*q                                                       //
#           W                a matrix of dimension n*(p*(q+1))                                               //
#           Y                a response vector of lentgh n                                                   //
#           X.test           a matrix of dimension ntest*q                                                   //
#           W.test           a matrix of dimension ntest*(p*(q+1))                                           //
#           Y.test           a response vector of lentgh ntest                                               //
#           alpha.initial    a scaler respesenting the initial value of intercept                            //
#           a0               a vector of length q, the initial value of a                                    //
#           b0               a vector of length p*(q+1), the initial value of b                              //                                         //
#           kappa1           the first tuning of prior method                                                //
#           kappa2           the second tuning of prior method                                               //
#           xi               the tuning of MCP                                                               //
#           epsilon          the tuning of Elastic net                                                       //
#           type             response type                                                                   //

# Output:                                                                                                    //
#           coef: a list of length three                                                                     //
#               -- alpha      a scaler respesenting the estimate of intercept                                //
#               -- a          the estimate vector of a                                                       //
#               -- b          the estimate vector of b                                                       //  
#           all.a           an array storing the estimate of a for all tuning combinations                   //                                       //
#           all.b           an array storing the estimate of b for all tuning combinations                   //                                       //
#           minmse          the minimum MSE value                                                            //
#           MSE             a matrix storing the mse value for all tuning combinations                       //
#           minEbic         the minimum MSE value                                                            //
#           Ebic            a matrix storing the EBIC value for all tuning combinations                      //
#           best.tuning      a vector of length two, representing the best tuning                            //
#-------------------------------------------------------------------------------------------------------------//


#-------------------------------------------------------------------------------------------------------------//
# Function: cv.sgMCP.differlam2
# Input:                                                                                                     //
#           X                a matrix of dimension n*q                                                       //
#           W                a matrix of dimension n*(p*(q+1))                                               //
#           Y                a response vector of lentgh n                                                   //
#           Yhat.prior       a predicted response vector                                                     //
#           X.test           a matrix of dimension ntest*q                                                   //
#           W.test           a matrix of dimension ntest*(p*(q+1))                                           //
#           Y.test           a response vector of lentgh ntest                                               //
#           alpha.initial    a scaler respesenting the initial value of intercept                            //
#           a0               a vector of length q, the initial value of a                                    //
#           b0               a vector of length p*(q+1), the initial value of b                              //                                         //
#           lambda1          a vector of the first tuning for psgMCP method                                  //
#           lambda2          a vector of the second tuning for psgMCP method                                 //
#           tau              a vector of the third for psgMPC method                                         //
#           xi               the tuning of MCP                                                               //
#           type             response type                                                                   //

# Output:                                                                                                    //
#           coef       a list of length three                                                                //
#               -- alpha      a scaler respesenting the estimate of intercept                                //
#               -- a          the estimate vector of a                                                       //
#               -- b          the estimate vector of b                                                       //
#           all.a            an array storing the estimate of a for all tuning combinations                  //                                       //
#           all.b            an array storing the estimate of b for all tuning combinations                  //                                       //
#           minmse           the minimum MSE value                                                           //
#           MSE              a matrix storing the mse value for all tuning combinations                      //
#           minEbic          the minimum MSE value                                                           //
#           Ebic             a matrix storing the EBIC value for all tuning combinations                     //
#           best.tuning      a vector of length two, representing the best tuning                            //
#------------------------------------------------------------------------------------------------------------//


#------------------------------------------------------------------------------------------------------------//
# Function: cv.sgMCP
# Input:                                                                                                     //
#           X                a matrix of dimension n*q                                                       //
#           W                a matrix of dimension n*(p*(q+1))                                               //
#           Y                a response vector of lentgh n                                                   //
#           X.test           a matrix of dimension ntest*q                                                   //
#           W.test           a matrix of dimension ntest*(p*(q+1))                                           //
#           Y.test           a response vector of lentgh ntest                                               //
#           alpha.initial    a scaler respesenting the initial value of intercept                            //
#           a0               a vector of length q, the initial value of a                                    //
#           b0               a vector of length p*(q+1), the initial value of b                              //                                         //
#           lambda1          a vector of the first tuning for psgMCP method                                  //
#           lambda2          a vector of the second tuning for psgMCP method                                 //
#           xi               the tuning of MCP                                                               //
#           type             response type                                                                   //

# Output:                                                                                                    //
#           coef       a list of length three                                                                //
#               -- alpha      a scaler respesenting the estimate of intercept                                //
#               -- a          the estimate vector of a                                                       //
#               -- b          the estimate vector of b                                                       //
#           all.a           an array storing the estimate of a for all tuning combinations                   //                                       //
#           all.b           an array storing the estimate of b for all tuning combinations                   //                                       //
#           minmse          the minimum MSE value                                                            //
#           MSE             a matrix storing the mse value for all tuning combinations                       //
#           minEbic         the minimum MSE value                                                            //
#           Ebic            a matrix storing the EBIC value for all tuning combinations                      //
#           best.tuning    a vector of length two, representing the best tuning                              //
#------------------------------------------------------------------------------------------------------------//

###############################################################################################################


###############################################################################################################

#--------------------------------               basic functions              ---------------------------------
taylor.exp  <-  function(X,W,Y,alpha0,a0,b0,type)
{
  if(type == "gaussian") {
    weight <- rep(1,n)
    Ystar <- Y
  }
  if(type=='poission') {
    eta <- X%*%a0+W%*%b0+alpha0
    mu <- exp(eta)
    weight <- mu
    Ystar <- eta+(Y-mu)/mu
  }
  return(list(weight=weight,Ystar=Ystar))
}

center.data <- function(X,W,Y,alpha0,a0,b0,type){
  taylor <- taylor.exp(X,W,Y,alpha0,a0,b0,type)
  weight <- taylor$weight
  Ystar <- taylor$Ystar
  mean.W <- apply(W,2,weighted.mean, w=weight)
  mean.Y <- weighted.mean(Ystar,weight)
  mean.X <- apply(X,2,weighted.mean, w=weight)
  meanlist <- list(mean.Y=mean.Y,mean.X=mean.X,mean.W=mean.W)
  Yols <- sqrt(weight)*(Ystar-mean.Y)
  weight.mat <- diag(as.vector(sqrt(weight)))
  Xols <- weight.mat%*%scale(X,center=mean.X,scale=FALSE)
  Wols <- weight.mat%*%scale(W,center=mean.W,scale=FALSE)
  return(list(Yols=Yols,Xols=Xols,Wols=Wols,meanlist=meanlist))
  
}

ARsigma <- function(p,rou){
  sigma <- diag(1,p)
  for(i in 1:p) {
    for(j in i:p) {
      sigma[i,j] <- rou^abs(i-j)
      sigma[j,i] <- sigma[i,j]
    }
  }
  return(sigma)
}

EBIC.fun <- function(X,W,Y,alpha,a,b,type) {
  r.EBIC <- 1-1/(2*log(p)/log(length(Y)))
  df <- length(which(c(a,b)!=0))
  eta.train <- alpha+X%*%a+W%*%b
  if(type=='gaussian') {
    error=Y-eta.train
    Ebic <- sum(error^2)+df*log(n)+2*df*r.EBIC*log(p)
    
  }
  if(type=='poisson') {
    logloss <- Y%*%eta.train-sum(exp(eta.train))
    Ebic <- -2*sum(logloss)+df*log(n)+2*df*r.EBIC*log(p)
  }
  
  return(Ebic=Ebic)
}

matW <- function(X,Z) {
  Xstar <- cbind(1,X)
  p <- ncol(Z)
  n <- nrow(X)
  W <- matrix(0,n,p*(q+1))
  
  for(j in 1:p) {
    ind <- ((q+1)*(j-1)+1):((q+1)*j)
    Zj <- Z[,j]
    Wj <- Zj*Xstar
    W[,ind] <- Wj
  }
  return(W = W)
}

generate.data <- function(alpha0,a,b,n,p,q,sigmaX,sigmaZ,sd.error,type) {
  X <- mvrnorm(n,rep(0,q),sigmaX)
  Z <- mvrnorm(n,rep(0,p),sigmaZ)
  W <- matW(X,Z)
  eta <- alpha0+X%*%a+W%*%b
  error <- rnorm(n,0,sd.error)
  if(type=='gaussian') {
    Y <- eta+error
  }
  
  if(type=='poisson') {
    X <- mvrnorm(n,rep(0.2,q),sigmaX)
    Z <- mvrnorm(n,rep(0.2,p),sigmaZ)
    W <- matW(X,Z)
    eta <- alpha0+X%*%a+W%*%b
    q.eta <- c(quantile(eta,0.02),quantile(eta,0.98))
    index <- which(eta>=q.eta[1]&eta<q.eta[2])
    eta <- eta[index]
    lambda <- exp(eta)
    Y <- rpois(length(lambda),lambda)
    X <- X[index,]
    Z <- Z[index,]
    W <- W[index,]
  }
  return(list(Y=Y,X=X,Z=Z,W=W)) 
  
}

matW.orth <- function(Wols) {
  W.orth <- matrix(0,nrow(Wols),p*(q+1))
  L.inverse.tran <- matrix(0,p*(q+1),p*(q+1))
  for(j in 1:p) {
    result <- matWj(Wols,j)
    ind <- ((q+1)*(j-1)+1):((q+1)*j)
    Wj.orth <- result$Wj.orth
    W.orth[,ind] <- Wj.orth
    L.inverse.tran[ind,ind]=result$L.inverse.tran
  }
  return(list(L.inverse.tran=L.inverse.tran,W.orth=W.orth))
}

orthogonized <- function(A) {
  L <- t(chol(A))
  L.tran.inverse <- solve(t(L))
  return(L.tran.inverse)
}

matWj <- function(Wols,j) {
  Wj <- matrix(0,nrow(Wols),q+1)
  indj <- ((q+1)*(j-1)+1):((q+1)*j)
  Wj <- Wols[,indj]
  n <- nrow(Wj)
  L.inverse.tran <- orthogonized((t(Wj)%*%Wj)/n)
  Wj.orth <- Wj%*%L.inverse.tran
  return(list(L.inverse.tran=L.inverse.tran,Wj.orth=Wj.orth))
}

vector.b <- function (beta,gamma) {
  mat.b.gamma <- cbind(beta,gamma)
  vec.b <- as.vector(t(mat.b.gamma))
  return (vec.b)
}


Sfun <- function(z,c) {
  znorm <- sqrt(sum(z^2))
  value <- max((1-c/znorm),0)*z
  return (value)
}

Sfun1<-function(z,c) {
  value<-max((1-c/abs(z)),0)*z
  return(value)
}

predict.sgMCP <- function(X.test,W.test,Y.test,alpha,abest,bbest,type) {
  n.test <- length(Y.test)
  eta.test <- as.vector(alpha)+X.test%*%abest+W.test%*%bbest
  if (type=='gaussian') {
    Y.pred <- eta.test
    error.test <- Y.test-eta.test
    mse <- sum(error.test^2)/n.test
    mape <- median(abs(error.test))
  }
  
  if (type=='poisson') {
    eta.true <- alpha+X.test%*%a+W.test%*%b
    lambda <- exp(eta.test)
    Y.pred <- lambda
    Y.pred <- ifelse(Y.pred>400,Y.test,Y.pred)
    Y.pred <- ifelse(abs(Y.pred-Y.test)>400,Y.test,Y.pred)
    mse <- sum(eta.true-eta.test)^2/n.test
    mape <- median(abs(Y.pred-Y.test))
  }
  
  return(list(mse=mse,Y.pred=Y.pred, MAPE=mape))
}

update.a <- function(Xols,Wols,Yols,alpha.initial,b0) {
  Ya <- Yols-Wols%*%b0
  fit <- lm(Ya~Xols)
  a0 <- coef(fit)[-1]
  alpha0 <- coef(fit)[1]
  return (list(a0=a0,alpha0=alpha0))
}
#---------------------------------------------------------------------------------------------------------


#-------------------------------               prior method           ------------------------------------
ridge.prior <- function(W,residual.prior,epsilon) {
  tw <- t(W)
  est <- (tw%*%residual.prior/n)/(1+epsilon)
  return(est)
}

GMCP.gamma <- function(Wj,bj,r,qgamma,kappa1,kappa2,xi) {
  muj <- (t(Wj)%*%r)/n+bj
  norm.bj <- sqrt(sum(bj^2))
  if(kappa1==0) {
    ghat <- 1
  } else {
    ghat <- ifelse (norm.bj==0,Inf,
                    1 + (1/norm.bj) * ifelse (norm.bj>xi*sqrt(qgamma)*kappa1,0,
                                              sqrt(qgamma)*kappa1-norm.bj/xi))
  }
  vj <- ifelse (abs(muj)>xi*kappa2*ghat,muj,
                sapply(muj,Sfun1,kappa2)/(1-1/(xi*ghat)))
  
  norm.vj <- sqrt(sum(vj^2))
  if(norm.vj > xi*sqrt(qgamma)*kappa1) {
    bjne <- vj
  } else {
    bjnew <- (xi/(xi-1))*Sfun(vj,sqrt(qgamma)*kappa1)
  }
  rnew <- r-Wj%*%(bjnew-bj)
  return (list(bj=bjnew,r=rnew))
}

GMCP.bj.prior.res <- function (Wj,bj,r,kappa1,kappa2,xi) {
  muj <- (t(Wj)%*%r)/n+bj
  norm.bj <- sqrt(sum(bj^2))
  if(kappa1==0){
    ghat <- 1
  } else {
    ghat <- ifelse(norm.bj==0,Inf,
                   1+(1/norm.bj)*ifelse(norm.bj>xi*sqrt(q+1)*kappa1,0,
                                        sqrt(q+1)*kappa1-norm.bj/xi))
  }
  vj <- muj
  vj[-1] <- ifelse(abs(muj[-1])>xi*kappa2*ghat,muj[-1],
                   sapply(muj[-1],Sfun1,kappa2)/(1-1/(xi*ghat)))
  norm.vj <- sqrt(sum(vj^2))
  if(norm.vj>xi*sqrt(q+1)*kappa1){
    bjnew <- vj
  } else{
    bjnew <- (xi/(xi-1))*Sfun(vj,sqrt(q+1)*kappa1)
  }
  rnew <- r-Wj%*%(bjnew-bj)
  return(list(bj=bjnew,r=rnew))
}

update.allb.prior <- function (W.orth,b0,r,kappa1,kappa2,xi,epsilon) {
  beta.ind <- (q+1)*(S_G-1)+1         #the index of elements in beta which are penalized using ridge
  jind <- S_EG[,1]
  kind <- S_EG[,2]
  gamma.ind <- (q+1)*(jind-1)+kind    #the index of elements in gamma which are penalized using ridge
  ridge.ind <- sort(c(beta.ind,gamma.ind))
  
  W.orth.ridge <- W.orth[,ridge.ind]
  b.ridge.old <- b0[ridge.ind]
  residual.prior <- r+W.orth.ridge%*%b.ridge.old
  b0[ridge.ind] <- coef(lm(residual.prior~W.orth.ridge))[-1]
  r <- residual.prior-W.orth.ridge%*%b0[ridge.ind]
  
  for(j in S_G){
    groupj <- setdiff((2:(q+1)),S_EG[which(S_EG[,1]==j),2])
    groupindex <- (q+1)*(j-1)+groupj
    Wj <- W.orth[,groupindex]
    bj <- b0[groupindex]
    qgamma <- length(groupindex)
    fit <- GMCP.gamma(Wj,bj,r,qgamma,kappa1,kappa2,xi)
    b0[groupindex] <- fit$bj
    r <- fit$r
  }
  
  bjind <- setdiff(seq(p),S_G)
  for(j in bjind) {
    jind <- ((j-1)*(q+1)+1):(j*(q+1))
    Wj <- W.orth[,jind]
    bj <- b0[jind]
    fit <- GMCP.bj.prior.res(Wj,bj,r,kappa1,kappa2,xi)
    b0[jind] <- fit$bj
    r <- fit$r
  }
  return(list(b=b0,r=r))
  
} 

update.para.prior  <-  function (Xols,Wols,Yols,alpha.initial,a0,b0,kappa1,kappa2,xi,epsilon) {
  a0 <- update.a (Xols,Wols,Yols,alpha.initial,b0) $a0
  orth <- matW.orth(Wols)
  r <- Yols-Xols%*%a0-Wols%*%b0
  W.orth <- orth$W.orth
  fit <- update.allb.prior(W.orth,b0,r,kappa1,kappa2,xi,epsilon)   #update the estimate of vector b
  b0 <- fit$b
  b0 <- orth$L.inverse.tran%*%b0
  return(list(a=a0,b=b0))
  
}

overall.estimate.prior <- function (X,W,Y,alpha.initial,a0,b0,type,kappa1,kappa2,xi,epsilon) {
  centdata <- center.data (X,W,Y,alpha.initial,a0,b0,type)
  Yols <- centdata$Yols
  Xols <- centdata$Xols
  Wols <- centdata$Wols
  mean.data <- centdata$meanlist
  
  est.new <- update.para.prior(Xols,Wols,Yols,alpha.initial,a0,b0,kappa1,kappa2,xi,epsilon)
  anew <- est.new$a
  bnew <- est.new$b
  alphanew <- mean.data$mean.Y-t(mean.data$mean.X)%*%anew-t(mean.data$mean.W)%*%bnew
  return(list(alphanew=alphanew,anew=anew,bnew=bnew))
}

sgMCP.prior <- function(X,W,Y,alpha.initial,a0,b0,type,kappa1,kappa2,xi,epsilon){
  s <- 0
  repeat{
    est.new <- overall.estimate.prior(X,W,Y,alpha.initial,a0,b0,type,kappa1,kappa2,xi,epsilon)
    aold <- a0
    a0 <- est.new$anew
    bold <- b0
    b0 <- est.new$bnew
    alphaold <- alpha.initial
    alpha.initial <- as.vector(est.new$alphanew)
    
    theta.new <- c(alpha.initial,a0,b0)
    theta.old <- c(alphaold,aold,bold)
    s <- s+1
    if(norm(theta.old-theta.new,'2') <= 0.001 || s>100 || max(abs(theta.new)) >= 3)
      break
    
  }
  return(list(alpha=alpha.initial,a=a0,b=b0))
}

cv.sgMCP.prior <- function (X,W,Y,X.test,W.test,Y.test,alpha.initial,a0,b0, type,kappa1,kappa2,xi,epsilon) {
  nkappa1 <- length(kappa1)
  nkappa2 <- length(kappa2)
  array.a <- array (0,dim=c(nkappa1,nkappa2,q)) 
  array.b <- array (0,dim=c(nkappa1,nkappa2,(q+1)*p))
  mat.alpha <- matrix(0,nrow=nkappa1,ncol=nkappa2)
  MSE <- matrix(0,nrow=nkappa1,ncol=nkappa2)
  Ebic <- matrix(0,nrow=nkappa1,ncol=nkappa2)
  
  for(i in 1:nkappa1)
  {
    for(j in 1:nkappa2)
    {
      fit <- sgMCP.prior(X,W,Y,alpha.initial,a0,b0,type,kappa1[i],kappa2[j],xi,epsilon)
      esta <- fit$a
      estb <- fit$b
      estalpha <- as.vector(fit$alpha)
      estlist <- fit
      a0 <- esta
      b0 <- estb
      alpha.initial <- estalpha
      
      prediction <- predict.sgMCP(X.test,W.test,Y.test,estalpha,esta,estb,type)
      MSE[i,j] <- prediction$mse
      Ebic[i,j] <- EBIC.fun(X,W,Y,estalpha,esta,estb,type)
      array.a[i,j,] <- esta
      array.b[i,j,] <- estb
      mat.alpha[i,j] <- estalpha
    }
  }
  
  matind <- which(Ebic==min(Ebic),arr.ind=TRUE)
  kappa1ind <- max(matind[,1])
  kappa2ind <- max(matind[which(matind[,1]==kappa1ind),2])
  best.tuning <- c(kappa1[kappa1ind],kappa2[kappa2ind])
  abest <- array.a[kappa1ind,kappa2ind,]
  bbest <- array.b[kappa1ind,kappa2ind,]
  alphabest <- mat.alpha[kappa1ind,kappa2ind]
  minmse <- min(MSE)
  minEbic <- min(Ebic)
  mse <- predict.sgMCP(X.test,W.test,Y.test,alphabest,abest,bbest,type)$mse
  
  return(list(coef=list(alpha=alphabest,a=abest,b=bbest),all.a=array.a,all.b=array.b,
              minmse=mse,MSE=MSE,minEbic=min(Ebic),Ebic=Ebic,best.tuning=best.tuning))  
}

#---------------------------------------------------------------------------------------------------------


#----------------------------                  psgMCP            ---------------------------------------
GMCP.bj.res <- function(Wj,bj,r,kappa1,kappa2,xi) {
  muj <- (t(Wj)%*%r)/n+bj
  norm.bj <- sqrt(sum(bj^2))
  if (kappa1==0) {
    ghat <- 1+epsilon.sgMCP
  } else {
    ghat <- ifelse(norm.bj==0,Inf,
                   1+epsilon.sgMCP+(1/norm.bj)*
                     ifelse(norm.bj>xi*sqrt(q+1)*kappa1,0, sqrt(q+1)*kappa1-norm.bj/xi))
  }
  vj <- muj
  vj[-1] <- ifelse (abs(muj[-1])>xi*kappa2*ghat,muj[-1],
                    sapply(muj[-1],Sfun1,kappa2)/(1-1/(xi*ghat)))
  norm.vj <- sqrt(sum(vj^2))
  if (norm.vj > xi*sqrt(q+1)*kappa1) {
    bjnew <- vj
  } else {
    bjnew <- (xi/((1+epsilon.sgMCP)*xi-1))*Sfun(vj,sqrt(q+1)*kappa1)
  }
  rnew <- r-Wj%*%(bjnew-bj)
  return(list(bj=bjnew,r=rnew))
}

update.allb <- function(W.orth,b0,r,lambda1,lambda2,xi) {
  for(j in 1:p) {
    jind <- ((j-1)*(q+1)+1):(j*(q+1))
    Wj <- W.orth[,jind]
    bj <- b0[jind]
    fit <- GMCP.bj.res(Wj,bj,r,lambda1,lambda2,xi)
    b0[jind] <- fit$bj
    r <- fit$r
  }
  return(list(b=b0,r=r))
}

update.para <- function(Xols,Wols,Yols,alpha,a0,b0,lambda1,lambda2,xi,epsilon) {
  a0 <- update.a(Xols,Wols,Yols,alpha,b0)$a0
  orth <- matW.orth(Wols)
  W.orth <- orth$W.orth
  r <- Yols-Xols%*%a0-Wols%*%b0
  fit <- update.allb(W.orth,b0,r,lambda1,lambda2,xi)###W的正交化放在内部
  b0 <- fit$b
  b0 <- orth$L.inverse.tran%*%b0
  
  return(list(a=a0,b=b0))
}

overall.estimate <- function (X,W,Y,alpha.initial,a0,b0,type,lambda1,lambda2,xi) {
  centdata <- center.data(X,W,Y,alpha.initial,a0,b0,type)
  Yols <- centdata$Yols
  Xols <- centdata$Xols
  Wols <- centdata$Wols
  mean.data <- centdata$meanlist
  est.new <- update.para(Xols,Wols,Yols,alpha.initial,a0,b0,lambda1,lambda2,xi)
  anew <- est.new$a
  bnew <- est.new$b
  alphanew <- mean.data$mean.Y-t(mean.data$mean.X)%*%anew-t(mean.data$mean.W)%*%bnew
  
  return(list(alphanew=alphanew,anew=anew,bnew=bnew))
}

sgMCP <- function(X,W,Y,alpha.initial,a0,b0,type,lambda1,lambda2,xi) {
  s <- 0
  repeat{
    est.new <- overall.estimate(X,W,Y,alpha.initial,a0,b0,type,lambda1,lambda2,xi)
    aold <- a0
    a0 <- est.new$anew
    bold <- b0
    b0 <- est.new$bnew
    alphaold <- alpha.initial
    alpha.initial <- as.vector(est.new$alphanew)
    
    
    theta.new <- c(alpha.initial,a0,b0)
    theta.old <- c(alphaold,aold,bold)
    s=s+1
    
    if (norm(theta.old-theta.new,'2') <=0.001 || s>100 || max(abs(theta.new))>=3)
      break
    
  }

  return(list(alpha=alpha.initial,a=a0,b=b0))
}

cv.sgMCP.differlam2 <- function(X,W,Y,Yhat.prior,X.test,W.test,Y.test,alpha.initial,a0,b0,
                                type,lambda1,lambda2,tau,xi) {
  ntau <- length(tau)
  nlambda1 <- length(lambda1)
  nlambda2 <- length(lambda2)
  
  array.a <- array(0,dim=c(ntau,nlambda1,nlambda2,q)) 
  array.b <- array(0,dim=c(ntau,nlambda1,nlambda2,(q+1)*p))
  array.alpha <- array(0,dim=c(ntau,nlambda1,nlambda2))
  MSE <- array(0,dim=c(ntau,nlambda1,nlambda2))
  Ebic <- array(0,dim=c(ntau,nlambda1,nlambda2))
  
  for(i in 1:ntau) {
    Ystar <- (Y+tau[i]*Yhat.prior)/(1+tau[i])
    for(j in 1:nlambda1) {
      for(k in 1:nlambda2) {
        fit <- sgMCP(X,W,Ystar,alpha.initial,a0,b0,type,lambda1[j],lambda2[k],xi)
        esta <- fit$a
        estb <- fit$b
        estalpha <- fit$alpha
        a0 <- esta
        b0 <- estb
        alpha.initial <- estalpha
        estlist <- fit
        
        prediction <- predict.sgMCP(X.test,W.test,Y.test,estalpha,esta,estb,type)
        MSE[i,j,k] <- prediction$mse
        Ebic[i,j,k] <- EBIC.fun(X,W,Y,estalpha,esta,estb,type)
        array.a[i,j,k,] <- esta
        array.b[i,j,k,] <- estb
        array.alpha[i,j,k] <- estalpha
      }
    }
  }
  matind <- which(Ebic==min(Ebic),arr.ind=TRUE)
  tauind <- max(matind[,1])
  lambda1ind <- max(matind[which(matind[,1]==tauind),2])
  lambda2ind <- max(matind[which(matind[,1]==tauind&matind[,2]==lambda1ind),3])
  
  best.tuning <- c(tau[tauind],lambda1[lambda1ind],lambda2[lambda2ind])
  abest <- array.a[tauind,lambda1ind,lambda2ind,]
  bbest <- array.b[tauind,lambda1ind,lambda2ind,]
  alphabest <- array.alpha[tauind,lambda1ind,lambda2ind]
  minmse <- min(MSE)
  
  return(list(coef=list(alpha=alphabest,a=abest,b=bbest),all.a=array.a,all.b=array.b,
              minmse=minmse,minEbic=min(Ebic),MSE=MSE,Ebic=Ebic,best.tuning=best.tuning))
  
}

#---------------------------------------------------------------------------------------------------------


#----------------------------                   sgMCP            ------------------------------------------
cv.sgMCP <- function(X,W,Y,X.test,W.test,Y.test,alpha.initial,a0,b0,
                     type,lambda1,lambda2,xi) {
  nlambda1 <- length(lambda1)
  nlambda2 <- length(lambda2)
  array.a <- array(0,dim=c(nlambda1,nlambda2,q)) 
  array.b <- array(0,dim=c(nlambda1,nlambda2,(q+1)*p))
  mat.alpha <- matrix(0,nrow=nlambda1,ncol=nlambda2)
  MSE=matrix(0,nrow=nlambda1,ncol=nlambda2)
  Ebic=matrix(0,nrow=nlambda1,ncol=nlambda2)
  
  for(i in 1:nlambda1)
  {
    for(j in 1:nlambda2)
    {
      
      fit <- sgMCP(X,W,Y,alpha.initial,a0,b0,type,lambda1[i],lambda2[j],xi)
      esta <- fit$a
      estb <- fit$b
      estalpha <- fit$alpha
      estlist <- fit###save as a list
      alpha.initial <- estalpha
      
      prediction <- predict.sgMCP(X.test,W.test,Y.test,estalpha,esta,estb,type)
      MSE[i,j] <- prediction$mse
      Ebic[i,j] <- EBIC.fun(X,W,Y,estalpha,esta,estb,type)
      array.a[i,j,] <- esta
      array.b[i,j,] <- estb
      mat.alpha[i,j] <- estalpha
    }
  }
  
  matind <- which(Ebic==min(Ebic),arr.ind=TRUE)
  lambda1ind <- max(matind[,1])
  lambda2ind <- max(matind[which(matind[,1]==lambda1ind),2])
  best.tuning <- c(lambda1[lambda1ind],lambda2[lambda2ind])
  abest <- array.a[lambda1ind,lambda2ind,]
  bbest <- array.b[lambda1ind,lambda2ind,]
  alphabest <- mat.alpha[lambda1ind,lambda2ind]
  minmse <- min(MSE)
  
  return(list(coef=list(alpha=alphabest,a=abest,b=bbest),all.a=array.a,all.b=array.b,minmse=minmse,
              MSE=MSE,minEbic=min(Ebic),Ebic=Ebic,best.tuning=best.tuning))
  
}

#---------------------------------------------------------------------------------------------------------


#############################                  the end               ######################################
