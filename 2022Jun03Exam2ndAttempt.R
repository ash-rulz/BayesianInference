####1(b)####
#Simulate 10000 draws from the posterior predictive distribution
nIter <- 10000
postTheta <- rgamma(nIter, 2326, 7)
postPredQ6 <- rpois(nIter, postTheta)

#2nd way - 1 posTheta, and 10k simulations of q6
postTheta1 <- rgamma(1, 2326, 7)
postPredQ61 <- rpois(nIter, postTheta)

#Plot a histogram of the draws
par(mfrow=c(2,1))
hist(postPredQ6)
hist(postPredQ61)
#Compute Pr(Q6 > 350|q1, ..., q5)
mean(postPredQ6 > 350)#0.1788

####1(c)####
#Use simulation to find the optimal a
aGrid <- seq(50, 1000, 1)
expUtil <- sapply(aGrid, function(a) mean(utility_func(a,postPredQ6)))
plot(aGrid, expUtil, type = 'l')
abline(v = aGrid[which.max(expUtil)])


####2(a)####
#Simulate 10000 draws from the joint posterior
set.seed(1)
mu_0 <- as.vector(rep(0,6))
Omega_0 <- (1/10^2)*diag(6)
v_0 <- 1
sigma2_0 <- 100^2
nIter <- 10000
postResult <- BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter)
postVar <- postResult$sigma2Sample
postBetaMat <- postResult$betaSample
#Posterior mean of all betas
postBetaMean <- apply(postBetaMat, 2, mean)
postBetaMean
#99% equal tail credible intervals for all parameters
apply(postBetaMat, 2, function(beta) quantile(beta, probs=c(0.005, 0.995)))
#99% probability that the b1 would lie between 5.5 and 12.7

####2(b)####
#Posterior mean of the sd
mean(sqrt(postVar))#40.00081
#Posterior median of the sd
median(sqrt(postVar))#39.57851

####2(c)####
x1 <- seq(min(X[,2]), max((X[,2])), 0.1)
x2 <- rep(27, length(x1))
const <- rep(1, length(x1))
newX <- cbind(const, x1, x1^2, x2, x2^2, x1*x2)
newMuMat <- postBetaMat %*% t(newX)
mu95Int <- apply(newMuMat, 2, 
                 function(mu) quantile(mu, probs=c(0.025, 0.975)))
plot(x1, mu95Int[1,], type ='l', ylim = c(0, 500))
lines(x1, mu95Int[2,], type ='l')


####2(d)####
x1 <- 50
x2 <- 25
const <- 1
newX <- cbind(const, x1, x1^2, x2, x2^2, x1*x2)
newMu <- as.vector(postBetaMat %*% t(newX))
newY <- newMu + rnorm(nIter, 0, sqrt(postVar))
hist(newY)

####2(e)####
tAct <- max(y)
newMuMat <- postBetaMat %*% t(X)
newY <- matrix(NA, nrow = nIter,ncol = length(y))
for (i in 1:length(y)) {
  newY[,i] <- newMuMat[,i] + rnorm(nIter, 0, sqrt(postVar))
}
tYRep <- apply(newY, 1, max)
mean(tYRep > tAct)#0.992


#####3(d)####
#Computes the log posterior distribution of theta
logPostTheta <- function(theta, x3Sum, n){
  logLike <- n*log(theta)-theta*x3Sum
  logPrior <- 2*log(theta)-4*theta
  return(logLike + logPrior)
}
x3Sum <- c(0.8, 1.1, 0.8, 0.9, 1)
x3Sum <- sum(x3Sum^3)
n <- 5
thetas <- seq(0.01, 2.5, 0.01)
#Get the unnormalized posterior of theta
postThetaUnnorm <- exp(logPostTheta(thetas, x3Sum, n))
#Normalize the postTheta distribution
normPostTheta <- postThetaUnnorm/(0.01*sum(postThetaUnnorm))

#Confirming that the area is 1
f <- approxfun(thetas, normPostTheta)
tot_area <- integrate(f, min(thetas), max(thetas))
print(paste("CDF:", tot_area$value))

plot(thetas, normPostTheta, type = 'l')

####3(e)####
initVal <- 0.1
OptimRes <- optim(initVal,
                  logPostTheta, lower=0.1,
                  gr=NULL,x3Sum, n,
                  method=c("L-BFGS-B"),
                  control=list(fnscale=-1),#Maximize log posterior
                  hessian=TRUE)#Returns the hessian
normApproxPost <- dnorm(thetas, mean = OptimRes$par,
                        sd = sqrt(solve(-OptimRes$hessian)))
lines(thetas, normApproxPost, col ='red')

