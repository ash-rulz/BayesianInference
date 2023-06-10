rm(list=ls())
set.seed(1)

####1(a)####
y1 <- c(2.32, 1.82, 2.40, 2.08,  2.13)
n <- length(y1)

####1(b)####
thetaPost <- rgamma(nIter, 2*n+1, 0.5+sum(y1))
nIter <- 10000
yPred <- rgamma(nIter, 2, thetaPost)
plot(density(yPred), type = 'l')
mean(yPred < 1.9)#0.706

####1(c)####
nWeeks <- c()
nIter <- 10000
for (i in 1:nIter) {
  #Draw theta for 30 weeks
  thetaPost <- rgamma(30, 11, 11.25)
  #No of weeks weights exceeds 2.4k KG
  nWeeks <- c(nWeeks,sum(sapply(thetaPost, 
                                function(theta) rgamma(1, 2, thetaPost)>2.4)))  
}
#hist(nWeeks)
mean(nWeeks)

####1(d)####
aGrid <- seq(2,10,0.01)#building cost
#aGrid <- c(1,2)
expLoss <- c()
for (a in aGrid) {
  breakW <- 0.9*log(a)#Weight at which escalator breaks
  #number of weeks out of the future 30 in which the escalator breaks
  nWeeks <- c()
  nIter <- 1000
  for (i in 1:nIter) {
    #Draw theta for 30 weeks
    thetaPost <- rgamma(30, 11, 11.25)
    #No of weeks weights exceeds breakW - the bridge breaks
    nWeeks <- c(nWeeks,sum(sapply(thetaPost, 
                                  function(theta) 
                                    rgamma(1, 2, thetaPost)>breakW)))  
  }
  loss <- a + nWeeks
  expLoss <- c(expLoss, mean(loss))
}
plot(aGrid, expLoss, type = 'l')
aGrid[which.min(expLoss)]



####2(a)####
mu_0 <- as.vector(rep(0,8))
Omega_0 <- (1/9)*diag(8)
v_0 <- 1
sigma2_0 <- 9
nIter <- 10000
X <- as.matrix(X)
postResult <- BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter)
postBetaMat <- postResult$betaSample
postVar <- postResult$sigma2Sample
quantile(postBetaMat[,2], probs = c(0.005, 0.995))
