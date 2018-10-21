#=====================================================================================================
# program: demo.R                                                                                     

# purpose: illustrate how to implement the proposed method with provided datasets

# steps: 
#    (1) source R files
#        source("functions.R")  
#    (2) generate datasets.  Demo datasets are provided for illustration with Simulation 1 under Scenario 1 incorporation prior set S_1.
#    (3) specify tunings
#    (4) set prior set. Demo prior set are provided for illustration with S_1
#    (5) implement the algorithm
#=====================================================================================================

# （1）source R files
library(MASS)
source("function.R")
start.time <- Sys.time()

#  (2) generate data
B <- 200#replication
n <- 200#sample size of training data
n.test <- 200#sample size of testing data
p <- 1000
q <- 4

sd.error <- 1
rouX <- 0.5
rouZ <- 0.3
sigmaX <- ARsigma(q,rouX)
sigmaZ <- ARsigma(p,rouZ)

type <- 'gaussian'
k.EBIC <- log(p)/log(n)
r.EBIC <- 1-1/(2*k.EBIC)

# coefficients and initial value
alpha0 <- 0
a <- runif(q,0.4,0.7)
beta <- c(runif(8,0.2,0.5),rep(0,p-8))
a0 <- rep(0,q)
b0 <- rep(0,(q+1)*p)
alpha.initial <- 0
gamma <- matrix(c(c(runif(3,0.2,0.5),0),
                  replicate(3,c(runif(1,0.2,0.5),rep(0,3))),
                  rep(0,p*q-16)),nrow=p,byrow=TRUE)
indexbeta <- seq(1,(q+1)*p,by=q+1)
b <- vector.b(beta,gamma)
vector.gamma <- b[-indexbeta]

#  (3) specify tunings
propor.tau <- seq(0,0.99,by=0.11)
tau <- propor.tau/(1-propor.tau)                                                                  # for psgMCP

lambda1.sgMCP <- 5^seq(-1.5,-0.56,by=0.02)                                                        # for sgMCP
lambda2.sgMCP <- 5^seq(-1.4,-0.56,by=0.02)                                                        # for sgMCP

lambda1 <- 5^seq(-1.5,-0.56,by=0.02)                                                              # for psgMCP
lambda2 <- 5^seq(-1.4,-0.56,by=0.02)                                                              # for psgMCP

kappa1 <- 5^seq(-1.5,-0.56,by=0.02)                                                               # for prior
kappa2 <- 5^seq(-1.4,-0.56,by=0.02)                                                               # for prior

xi <- 6                                                                                           # for MCP
epsilon <- 0.05                                                                                   # for ridge in prior method
epsilon.sgMCP <- 0.05                                                                             # for ridge in sgMCP

# (4) prior set S_1:  7T+2F+4T+2F
S_G0 <- c(5:8,20)                         
S_EG <- cbind(c(2:4,22,22,2),c(rep(2,4),3,3))                                                     # for E-G interactions
S_H <- unique(S_EG[,1])              
S_G <- sort(union(S_G0,S_H))                                                                      # for G variables


#  (5) implement the algprithm   
TP.main <- matrix(0,ncol=3,nrow=B)                                                                # for main effects
colnames(TP.main) <- c('TP.sgMCP.main','TP.prior.main','TP.psgMCP.main')
FP.main <- matrix(0,ncol=3,nrow=B)                                                                # for main effects
colnames(FP.main) <- c('FP.sgMCP.main','FP.prior.main','FP.psgMCP.main')

MSE.whole <- matrix(0,ncol=3,nrow=B)                                                               # for whole model
colnames(MSE.whole) <- c('MSE.sgMCP.whole','MSE.prior.whole','MSE.psgMCP.whole')
MSE.main <- matrix(0,ncol=3,nrow=B)                                                                # for main effects
colnames(MSE.main) <- c('MSE.sgMCP.main','MSE.prior.main','MSE.psgMCP.main')
Bias.beta <- matrix(0,ncol=3,nrow=B)                                                                # for main effects
colnames(Bias.beta) <- c('Bias.sgMCP.main','Bias.prior.main','Bias.psgMCP.main')


TP.interaction <- matrix(0,ncol=3,nrow=B)                                                            # for interaction effects
colnames(TP.interaction) <- c('TP.sgMCP.interaction','TP.prior.interaction',
                              'TP.psgMCP.interaction')
FP.interaction <- matrix(0,ncol=3,nrow=B)                                                            # for interaction effects
colnames(FP.interaction) <- c('FP.sgMCP.interaction','FP.prior.interaction','FP.psgMCP.interaction')

MSE.inter <- matrix(0,ncol=3,nrow=B)                                                                 # for interaction effects
colnames(MSE.inter) <- c('MSE.sgMCP.inter','MSE.prior.inter','MSE.psgMCP.inter')  
Bias.gamma <- matrix(0,ncol=3,nrow=B)                                                                # for interaction effects
colnames(Bias.gamma) <- c('Bias.sgMCP.interaction','Bias.prior.interaction',
                          'Bias.psgMCP.interaction')


tuning.sgMCP <- matrix(0,nrow=B,ncol=2)                                                               # for sgMCP
colnames(tuning.sgMCP) <- c("lambda1.sgMCP",'lambda2.sgMCP')
tuning.prior <- matrix(0,nrow=B,ncol=2)                                                               # for prior
colnames(tuning.prior) <- c("kappa1.prior",'kappa2.prior')
tuning.best <- matrix(0,nrow=B,ncol=4)                                                                # for psgMCP
colnames(tuning.best) <- c("tau.psgMCP",'lambda1.psgMCP','lambda2.psgMCP','propor.tau')


#---------------------         begin B replications            ---------------------
for (i in 1:B) {
  
  #---------------------------    1. data     ---------------------------------------
  #---------------------------      1.1 training data (X,W,Y)     -------------------
  simulation1.data <- generate.data(alpha0=0,a,b,n,p,q
                                  ,sigmaX,sigmaZ,sd.error,type=type)
  X <- simulation1.data$X
  Z <- simulation1.data$Z
  W <- simulation1.data$W
  Y <- simulation1.data$Y

  #-----------------------      1.2 testing data (X.test,W.test,Y,test)     ----------------
  simulation1.testdata <- generate.data(alpha0=0,a,b,n.test,p,q
                                      ,sigmaX,sigmaZ,sd.error,type=type)
  X.test <- simulation1.testdata$X
  Z.test <- simulation1.testdata$Z
  W.test <- simulation1.testdata$W
  Y.test <- simulation1.testdata$Y
  

  #---------------------                2. prior method                -----------------------
  fit.prior <- cv.sgMCP.prior(X,W,Y,X.test,W.test,Y.test,alpha.initial,a0,b0,
                            type,kappa1,kappa2,xi,epsilon)
  coef.prior <- fit.prior$coef
  a.prior <- coef.prior$a
  b.prior <- coef.prior$b
  alpha.prior <- coef.prior$alpha
  
  index.b.prior <- which(b.prior!=0)
  W.refit <- W[,index.b.prior]
  refit.prior <- lm(Y~X+W.refit)
  coef.refit <- coef(refit.prior)[-1]
  a.refit.prior <- coef.refit[1:q]
  b.refit.prior <- coef.refit[-(1:q)]
  b.prior[index.b.prior] <- b.refit.prior
  
  #--------------------                 3. estimate of Yhat.prior         -----------------------
  priorfit <- predict.sgMCP(X,W,Y,alpha.prior,a.refit.prior,b.prior,type)
  Yhat.prior <- priorfit$Y.pred

  #--------------------                 4.psgMCP                        --------------------------
  
  fit.final <-  cv.sgMCP.differlam2(X,W,Y,Yhat.prior,X.test,W.test,Y.test,alpha.initial,a0,b0,
                                  type,lambda1,lambda2,tau,xi)
  
  coef.final <- fit.final$coef
  a.final <- coef.final$a
  b.final <- coef.final$b
  alpha.final <- coef.final$alpha
  
  
  #------------------                    5.sgMCP                        --------------------------
  fit.sgMCP <- cv.sgMCP(X,W,Y,X.test,W.test,Y.test,alpha.initial,a0,b0,
                      type,lambda1.sgMCP,lambda2.sgMCP,xi)
  coef.sgMCP <- fit.sgMCP$coef
  a.sgMCP <- coef.sgMCP$a
  b.sgMCP <- coef.sgMCP$b
  alpha.sgMCP <- coef.sgMCP$alpha
   

  #-------------------                  6. Summary of results           --------------------------
  
  #-------------------                  TP and FP for main effects      ---------------------------
  beta.sgMCP <- b.sgMCP[indexbeta]
  TP.beta.sgMCP <- sum((beta!=0)&beta.sgMCP!=0)
  FP.beta.sgMCP <- sum((beta==0)&beta.sgMCP!=0)
  TP.main[i,1] <- TP.beta.sgMCP
  FP.main[i,1] <- FP.beta.sgMCP
  
  beta.prior <- b.prior[indexbeta]
  TP.beta.prior <- sum((beta!=0)&beta.prior!=0)
  FP.beta.prior <- sum((beta==0)&beta.prior!=0)
  TP.main[i,2] <- TP.beta.prior
  FP.main[i,2] <- FP.beta.prior
  
  beta.final <- b.final[indexbeta]
  TP.beta <- sum((beta!=0)&beta.final!=0)
  FP.beta <- sum((beta==0)&beta.final!=0)
  TP.main[i,3] <- TP.beta
  FP.main[i,3] <- FP.beta
  
 
  #--------------------             TP and FP for interactions     --------------------------
  vector.gamma.sgMCP <- b.sgMCP[-indexbeta]
  TP.gamma.sgMCP <- sum((vector.gamma!=0)&vector.gamma.sgMCP!=0)
  FP.gamma.sgMCP <- sum((vector.gamma==0)&vector.gamma.sgMCP!=0)
  TP.interaction[i,1] <- TP.gamma.sgMCP
  FP.interaction[i,1] <- FP.gamma.sgMCP
  
  vector.gamma.prior <- b.prior[-indexbeta]
  TP.gamma.prior <- sum((vector.gamma!=0)&vector.gamma.prior!=0)
  FP.gamma.prior <- sum((vector.gamma==0)&vector.gamma.prior!=0)
  TP.interaction[i,2] <- TP.gamma.prior
  FP.interaction[i,2] <- FP.gamma.prior
  
  vector.gamma.final <- b.final[-indexbeta]
  TP.gamma <- sum((vector.gamma!=0)&vector.gamma.final!=0)
  FP.gamma <- sum((vector.gamma==0)&vector.gamma.final!=0)
  TP.interaction[i,3] <- TP.gamma
  FP.interaction[i,3] <- FP.gamma

  #-------------------                         PMSE             ------------------------------
  ## for sgmcp
  index.b.sgMCP <- which(b.sgMCP!=0)
  W.refit <- W[,index.b.sgMCP]
  refit.sgMCP <- lm(Y~X+W.refit)
  coef.refit <- coef(refit.sgMCP)[-1]
  a.refit.sgMCP <- coef.refit[1:q]
  b.refit.sgMCP <- coef.refit[-(1:q)]
  b.sgMCP[index.b.sgMCP] <- b.refit.sgMCP
  sgMCPfit <- predict.sgMCP(X.test,W.test,Y.test,alpha.sgMCP,a.refit.sgMCP,b.sgMCP,type)
  
  MSE.sgMCP <- sgMCPfit$mse
  MSE.whole[i,1] <- MSE.sgMCP
  MSE.main[i,1] <- sum((W.test[,indexbeta]%*%(b.sgMCP[indexbeta]-beta))^2)/n.test# mse of main part
  MSE.inter[i,1] <- sum((W.test[,-indexbeta]%*%(b.sgMCP[-indexbeta]-vector.gamma))^2)/n.test# mse of main part
  
  
  # for prior method
  priorfit <- predict.sgMCP(X.test,W.test,Y.test,alpha.prior,a.refit.prior,b.prior,type)
  MSE.prior <- priorfit$mse
  MSE.whole[i,2] <- MSE.prior
  MSE.main[i,2] <- sum((W.test[,indexbeta]%*%(b.prior[indexbeta]-beta))^2)/n.test# mse of main part
  MSE.inter[i,2] <- sum((W.test[,-indexbeta]%*%(b.prior[-indexbeta]-vector.gamma))^2)/n.test# mse of main part
  
  
  #for psgMCP
  index.b.final <- which(b.final!=0)
  W.refit <- W[,index.b.final]
  refit.final <- lm(Y~X+W.refit)
  coef.refit <- coef(refit.final)[-1]
  a.refit.final <- coef.refit[1:q]
  b.refit.final <- coef.refit[-(1:q)]
  b.final[index.b.final] <- b.refit.final
  
  finalfit <- predict.sgMCP(X.test,W.test,Y.test,alpha.final,a.refit.final,b.final,type)
  MSE.final <- finalfit$mse
  MSE.whole[i,3] <- MSE.final
  MSE.main[i,3] <- sum((W.test[,indexbeta]%*%(b.final[indexbeta]-beta))^2)/n.test
  MSE.inter[i,3] <- sum((W.test[,-indexbeta]%*%(b.final[-indexbeta]-vector.gamma))^2)/n.test
  
  
  #--------------------------                            Bias                 ------------------------------
  Bias.beta[i,3] <- sum(abs(beta-beta.final))/sum(abs(beta))                                       # for psgMCP
  Bias.beta[i,2] <- sum(abs(beta-beta.prior))/sum(abs(beta))                                       # for prior
  Bias.beta[i,1] <- sum(abs(beta-beta.sgMCP))/sum(abs(beta))                                       # for sgMCP
  
  Bias.gamma[i,3] <- sum(abs(vector.gamma-vector.gamma.final))/sum(abs(vector.gamma))              # for psgMCP
  Bias.gamma[i,2] <- sum(abs(vector.gamma-vector.gamma.prior))/sum(abs(vector.gamma))              # for prior
  Bias.gamma[i,1] <- sum(abs(vector.gamma-vector.gamma.sgMCP))/sum(abs(vector.gamma))              # for sgMCP
  
  
  #------------------------                          tuning                 -------------------------------
  tuning.best[i,1:3] <- fit.final$best.tuning                                                      # for psgMCP
  tuning.best[i,4] <-  tuning.best[i,1]/(1+ tuning.best[i,1])
  tuning.prior[i,] <- fit.prior$best.tuning                                                        # for prior
  tuning.sgMCP[i,] <- fit.sgMCP$best.tuning                                                        # for sgMCP
  
  message('This is the ',i,'th ','replication Sm1.1-S1')
  allresult <- cbind(TP.main,FP.main,TP.interaction,FP.interaction,MSE.whole,MSE.main,MSE.inter,
                   Bias.beta,Bias.gamma,tuning.best,tuning.sgMCP,tuning.prior)
  print(allresult)
  write.csv(allresult,
            file='Sm1.1-S1.csv')
}

allresult <- cbind(TP.main,FP.main,TP.interaction,FP.interaction,MSE.whole,MSE.main,MSE.inter,
                 Bias.beta,Bias.gamma,tuning.best,tuning.sgMCP,tuning.prior)
all.median <- round(apply(allresult,2,median),2)
all.mean <- round(apply(allresult,2,mean),2)
all.sd <- round(apply(allresult,2,sd),2)
summary.result <- rbind(all.median,all.mean,all.sd)
row.names(summary.result) <- c('median','mean','sd')
write.csv(allresult, file = 'Sm1.1-S1.csv')
write.csv(summary.result, file = 'summary for Sm1.1-S1.csv')


end.time <- Sys.time()-start.time
end.time
