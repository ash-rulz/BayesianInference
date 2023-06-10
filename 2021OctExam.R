####1(d)####
logPostFn <- function(theta, x, n){
  logLik <- x * log(theta) - n * theta
  logPrior <- 2*log(theta) - 0.5*theta
  return(logLik + logPrior)
}
thetas <- seq(3, 8, 0.01)
x <- 75
n <- 15
logPostUnNorm <- sapply(thetas, function(theta) exp(logPostFn(theta, x, n)))
logPostNorm <- logPostUnNorm / (0.01 * sum(logPostUnNorm))
plot(thetas, logPostNorm, type = 'l', 
     ylab = 'P(theta|X)')

####1(e)####
initVal <- 3.2
OptimRes <- optim(initVal,
                  logPostFn,
                  gr=NULL,x, n, 
                  lower=3,
                  method=c("L-BFGS-B"),
                  control=list(fnscale=-1),#Maximize log posterior
                  hessian=TRUE)#
normPostSamples <- sapply(thetas, 
                          function(theta) rnorm(1, OptimRes$par,
                                            sqrt(solve(-OptimRes$hessian))))
normPostSamples <- rnorm(10000, OptimRes$par,
      sqrt(solve(-OptimRes$hessian)))
lines(density(normPostSamples), col = 'red')

####1(f)####
#Draw from posterior Gamma(78,15.5)
tRep <- c()
for (i in 1:10000) {
  thetas <- rgamma(15, 78, 15.5)
  tRep <- c(xMaxVector, 
                  max(sapply(thetas, function(theta) rpois(1, theta))))  
}
mean(tRep >= 14)

####2(a)####
nIter <- 10000
mu_0 <- as.vector(rep(0,3))
Sigma_0 <- 16*diag(3)
betaPostMatrix <- BayesLogitReg(y, X, mu_0, Sigma_0, nIter)$betaSample
#Compute the 90% equal tail credible interval for β1 and interpret it
beta1Interval <- quantile(betaPostMatrix[,2], probs = c(0.05, 0.95))

####2(b)####
#Compute the posterior probability that β2 > 0
beta2Prob <- mean(betaPostMatrix[,3] > 0)#0.8882


####2(c)####
#Compute the joint posterior probability that both β1 > 0 and β2 > 0
beta13Prob <- mean(betaPostMatrix[,3] > 0 & 
                    betaPostMatrix[,2] > 0)#0.8747

####2(d)####
probVect <- sapply(betaPostMatrix[,1], 
                   function(beta0) exp(beta0)/(1+exp(beta0)))
mean(probVect>0.5)#0.0136

####2(e)####
x1 <- seq(min(X[,2]), max(X[,2]), 0.01)
const <- rep(1, length(x1))
x2 <- rep(1, length(x1))
newX <- cbind(const, x1, x2)
linPred <- betaPostMatrix %*% t(newX)
pKMat <- exp(linPred)/(1+exp(linPred))
pk95Int <- apply(pKMat, 2, 
                 function(x) quantile(x, probs=c(0.025, 0.975)))
plot(x1, pk95Int[1,], type = 'l', ylim = c(0,1))
lines(x1, pk95Int[2,], type = 'l')


####3(b)####
#Simulate mu_n from posterior ∼ N(92, 4)
mu_n <- rnorm(10000, 92, 2)
postPredVect <- sapply(mu_n, function(mu) rnorm(1, mu, sqrt(50)))
hist(postPredVect)

####3(c)####

adCost <- seq(0, 20, 0.01)
util <- c()
rev <- c()
for (i in 1:length(adCost)) {
  #Draw 10000 mu from posterior
  muN <- rnorm(10000, 92, 2)
  #For each muN calc the utility
  util <- c(util, 
            mean(sapply(muN, 
                        function(mu)(60+sqrt(adCost[i])*log(mu)-adCost[i]))))
}
par(mfrow=c(1,1))
plot(adCost, util, type = 'l', col = 'green')