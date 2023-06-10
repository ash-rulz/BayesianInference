####1(a)####
#Compute the posterior probability that θA > 0.4
#Draw from posterior Beta(54,86)
postThetaA <- rbeta(10000, 54, 86)
mean(postThetaA > 0.4)#0.3627

#plot the posterior distribution of 1 − θA
hist((1-postThetaA), main = 'Histogram of probability of choosing non A brands')

####1(b)####
#Compute a 95% equal tail credible interval for the ratio 1−θA/θA and interpret it
#(probability of choosing brand B & C):(the probability of choosing brand A)
postNonAToARatio <- (1-postThetaA)/postThetaA
hist(postNonAToARatio)
quantile(postNonAToARatio, probs = c(0.025, 0.975))#1.141642 2.270124
#The 95% equi tailed credible interval of the customer preferring
#brand B and C over A is 1.141642 2.270124

####1(c)####
#Marginal likelihood of the model
marginalLhood <- beta(54,86) / beta(16,24)
marginalLhood

####1(d)####
priorAlpha <- c(20,20,20)
x <- c(38,27,35)
thetaA <- rgamma(10000, shape=priorAlpha[1]+y[1],rate= 1)
thetaB <- rgamma(10000, shape=priorAlpha[2]+y[2],rate= 1)
thetaC <- rgamma(10000, shape=priorAlpha[3]+y[3],rate= 1)
thetaMat <- cbind(thetaA, thetaB, thetaC)
for (i in 1:nrow(thetaMat)) {
  total <- sum(thetaMat[i,])
  thetaMat[i,1] <- thetaMat[i,1]/total
  thetaMat[i,2] <- thetaMat[i,2]/total
  thetaMat[i,3] <- thetaMat[i,3]/total
}
mean(as.vector(thetaMat[,1]) > as.vector(thetaMat[,3]))

###2(d)####
#Computes the log posterior distribution of θ
logPosteriorFn <- function(theta, x2Sum, n){
  return((n*log(theta)) - 
           (theta*(x2Sum + 0.5)))
}
thetas <- seq(0.1, 10, 0.01)
x2Sum <- 2.8
n <- 13
posPdf <- exp(logPosteriorFn(thetas, x2Sum, n))
normPosPdf <- posPdf/(0.01*sum(posPdf))
plot(thetas, normPosPdf, type = 'l', ylab = 'Posterior Distribution')


###2(e)####
initVal <- 1
OptimRes <- optim(initVal,
                  logPosteriorFn,
                  gr=NULL,x2Sum, n,
                  lower=0.1, 
                  method=c("L-BFGS-B"),
                  control=list(fnscale=-1),#Maximize log posterior
                  hessian=TRUE)#Returns the hessian
postSamples <- rnorm(10000, mean = OptimRes$par, 
                     sd = sqrt(solve(-OptimRes$hessian)))
lines(density(postSamples), col = 'red')


####3(a)####
#Simulate 10000 draws from the joint posterior
mu_0 <- rep(0,7)
Omega_0 <- 1/25*diag(7)
v_0 <- 1
sigma2_0 <- 4
nIter <- 10000
postDraws <- BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter)
postBetas <- postDraws$betaSample
postVar <- postDraws$sigma2Sample
#Compute the posterior mean
apply(postBetas, 2, mean)
#1.30890328  0.69975799  0.15738372  0.42682638 -0.16264437  0.07898822 -0.24059971

#95% equal tail credible intervals for all parameters
apply(postBetas, 2, function(x) quantile(x, probs=c(0.025, 0.975)))

####3(b)####
#Compute the posterior median of the standard deviation σ
median(sqrt(postVar))#0.6395435

#Investigate the effect of x1*x3 & x1*x4
plot(X[,6], X[,7])
postY <- postBetas %*% t(X)
for (i in 1:ncol(postY)) {
  postY[,i] <- postY[,i] + 
    sapply(postVar, function(x) rnorm(1, 0, sqrt(x)))
}
#Ummbi

####3(d)####
newBetas <- postBetas[,1:3]
x1 <- seq(min(X[,2]), max(X[,2]), by = 0.01)
const <- rep(1, length(x1))
x2 <- rep(0.5, length(x1))
newX <- cbind(const, x1, x2)
newMu <- newBetas %*% t(newX)
mu90Interval <- apply(newMu, 
                      2,
                      function(x) quantile(x, probs=c(0.05, 0.95)))
plot(x1,mu90Interval[1,],type = 'l')
lines(x1,mu90Interval[2,],type = 'l')


####3(e)####
x1 <- 0.4
x2 <- 1
x3 <- 1
x4 <- 0
newX <- matrix(c(1, x1, x2, x1 * x3), ncol = 1)
newBetas <- postBetas[,1:3]
newBetas <- cbind(newBetas,
                  postBetas[,6])
newMu <- as.vector(newBetas %*% newX)
newY <- newMu + sapply(postVar, function(x) rnorm(1, 0, sqrt(x)))
hist(newY)
