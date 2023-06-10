####2(b)####
13/(2.8+.5)

####2(d)####
logPostFn <- function(theta, xSumSq, n){
  logLik <- n * log(theta) -
    theta * xSumSq
  logPrior <- -0.5*theta
  logPost <- logLik + logPrior
  return(logPost)
}
thetaGrid <- seq(0.5, 9, 0.01)
xSumSq <- 2.8
n <- 13
#Unnormalized theta posterior
thetaPostUn <- exp(logPostFn(thetaGrid, xSumSq, n))
thetaPostNorm <- thetaPostUn/(0.01 * sum(thetaPostUn))

#Confirming that the area is 1
f <- approxfun(thetaGrid, thetaPostNorm)
tot_area <- integrate(f, min(thetaGrid), max(thetaGrid))
print(paste("CDF:", tot_area$value))

plot(thetaGrid, thetaPostNorm, type = 'l')

####2(e)####
initVal <- 0.5
OptimRes <- optim(initVal,
                  logPostFn,lower=0.1, 
                  gr=NULL, xSumSq = xSumSq, n = n,
                  method=c("L-BFGS-B"),
                  control=list(fnscale=-1),#Maximize log posterior
                  hessian=TRUE)#Returns the hessian
thetaApprox <- dnorm(thetaGrid, OptimRes$par, sqrt(solve(-OptimRes$hessian)))
lines(thetaGrid, thetaApprox, col = 'red')


####1(a)####
nIter <- 10000
thetaAPost <- rbeta(nIter, 54, 86)
mean(thetaAPost > 0.4)#0.359

thetaBCPost <- 1-thetaAPost
hist(thetaBCPost)

####1(b)####
#Odds of not choosing brand A
thetaOdds <- (1-thetaAPost)/thetaAPost
hist(thetaOdds)
quantile(thetaOdds, probs = c(0.025, 0.975))
#There is relatively high odds of not choosing brand A

####1(c)####
beta(54,86)/beta(16,24)


####

####3(a)####
mu_0 <- rep(0, 7)
Omega_0 <- (1/5^2)*diag(7)
v_0 <- 1
sigma2_0 <- 4
nIter <- 10000
postResult <- BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter)
postBetaMat <- postResult$betaSample
postVar <- postResult$sigma2Sample

apply(postBetaMat, 2, mean)
apply(postBetaMat, 2, 
      function(beta) quantile(beta, probs=c(0.025, 0.975)))

####3(b)####
median(sqrt(postVar))

####3(c)####


####3(d)####
x1 <- seq(min(X[,2]), max(X[,2]), 0.01)
const <- rep(1, length(x1))
x2 <- rep(0.5, length(x1))
x3 <- rep(0,length(x1))
x4 <- rep(0,length(x1))
newX <- cbind(const, x1, x2, x3, x4, x1 * x3, x1 * x4)
newMu <- postBetaMat %*% t(newX)
newMuInt <- apply(newMu, 2, function(x) quantile(x, probs=c(0.05, 0.95)))
plot(x1, newMuInt[1,], type = 'l')
lines(x1, newMuInt[2,], type = 'l')

####3(e)####
x1 <- 0.4
x2 <- x3 <- const <- 1
x4 <- 0
newX <- c(const, x1, x2, x3, x4, x1 * x3, x1 * x4)
newMu <- postBetaMat %*% newX
newY <- sapply(1:nIter, 
       function(x) rnorm(1, mean = newMu[x], sd = sqrt(postVar[x])))
hist(newY)
