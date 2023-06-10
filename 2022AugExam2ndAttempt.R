####1(b)####
x <- c(0.7, 1.1, 0.9, 1.5)
prod(x^4)
n <- 4
(gamma(4+3*4)*2^4)/
  (gamma(4)*(2+sum(1/x))^(4+3*n)*2^n*prod(x^4))

thetaHat <- (3*n+3)/(2+sum(1/x))
term1 <- 3*n*log(thetaHat) - thetaHat*sum(1/x) - n*log(2) -
  4*sum(log(x))
a <- 4
b <- 2
therm2 <- a*log(b) + 3*log(thetaHat) - 2*thetaHat - 
  log(gamma(a))
term3 <- log(thetaHat^2/(3*n+3))/2
exp(term1+therm2+term3+(log(2*pi)/2))

####1(d)####
logPostFn <- function(theta, x){
  n <- length(x)
  logLik <- 3*n*log(theta) - theta*sum(1/x)
  logPrior <- 3*log(theta) - 2*theta
  return(logLik + logPrior)
}
x <- c(0.7, 1.1, 0.9, 1.5)
thetaGrid <- seq(0.5, 6, 0.01)
postThetaUnnorm <- exp(logPostFn(thetaGrid, x))
postThetaNorm <- postThetaUnnorm/(0.01*sum(postThetaUnnorm))

f <- approxfun(thetaGrid, postThetaNorm)
tot_area <- integrate(f, min(thetaGrid), max(thetaGrid))
print(paste("CDF:", tot_area$value))

plot(thetaGrid, postThetaNorm, type = 'l', 
     xlab = 'Theta'
     ,ylab = 'Normalized Posterior Density')

####1(e)####
initVal <- 0.1
OptimRes <- optim(initVal,
                  logPostFn, lower=0.1, 
                  gr=NULL, x,
                  method=c("L-BFGS-B"),
                  control=list(fnscale=-1),#Maximize log posterior
                  hessian=TRUE)#Returns the hessian
normApproxThetaDensity <- dnorm(thetaGrid,
                                mean = OptimRes$par,
                                sd = sqrt(solve(-OptimRes$hessian)))
lines(thetaGrid, normApproxThetaDensity, col = 'red')


####2(a)####
#Simulate 20000 draws from the joint posterior
mu_0 <- as.vector(rep(0, 3))
Sigma_0 <- 100*diag(3) 
nIter <- 20000
postBetaMat <- (BayesLogitReg(y, X, mu_0, Sigma_0, nIter))$betaSample
#Compute a 95% equal tail credible interval for β1
quantile(postBetaMat[,2], probs = c(0.025, 0.975))
#     2.5%      97.5% 
#0.01351641 0.18288872 

####2(b)####
#Compute the joint posterior probability that both β1 > 0 and β2 > 0.
mean((postBetaMat[,2] > 0) & (postBetaMat[,3] > 0))#0.9204

####2(c)####
x1 <- 5
x2 <- 1
newX <- matrix(c(1, x1, x2), nrow = 3)
linPred <- postBetaMat %*% newX
p <- exp(linPred)/(1+exp(linPred))
#Odds of not repairing the bridge
pOdds <- (1-p)/p
hist(pOdds)


####2(d)####
x1 <- seq(min(X[,2]), max(X[,2]), 0.1)
x2 <- rep(0, length(x1))
const <- rep(1, length(x1))
newX <- cbind(const, x1, x2)
linPred <- postBetaMat %*% t(newX)
pMat <- exp(linPred)/(1+exp(linPred))
pInt <- apply(pMat, 2,
      function(x) quantile(x, probs = c(0.025, 0.975)))
plot(x1, pInt[1,], type = 'l', ylim = c(0,1))
lines(x1, pInt[2,], type = 'l')


####2(e)####
x1 <- 40
x2 <- 1
newX <- matrix(c(1, x1, x2), nrow = 3)
linPred <- postBetaMat %*% newX
p <- exp(linPred)/(1+exp(linPred))
hist(p)
mean(p>0.5)#0.1621

beta(28,22)/beta(39,22)
