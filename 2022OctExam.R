####1(d)####
#Function for log posterior distribution of θ as the sum
#of the log likelihood and the log prior distribution of θ
logPostFn <- function(theta){
  x <- c(13, 8, 11, 7)
  n <- length(x)
  N <- 20
  logLikelihood <- sum(x) * log(theta) + 
    (n*N - sum(x)) * log(1-theta)
  logPrior <- log(theta) + 2*log(1-theta)
  return(exp(logLikelihood + logPrior))
}
thetas <- seq(0,1,0.001)

#Calculating the normalizing constant for given unknown distribution
normConst <- integrate(logPostFn, lower = min(thetas), 
                        upper = max(thetas))
normPostDistr <- logPostFn(thetas)/normConst$value
plot(thetas, normPostDistr, type = 'l', xlab = 'theta'
     ,ylab = 'Normalized Posterior Density')

#Confirming that the area is 1
f <- approxfun(thetas, normPostDistr)
tot_area <- integrate(f, min(thetas), max(thetas))
print(paste("CDF:", tot_area$value))


####1(e)####
logPostFn <- function(theta){
  x <- c(13, 8, 11, 7)
  n <- length(x)
  N <- 20
  logLikelihood <- sum(x) * log(theta) + 
    (n*N - sum(x)) * log(1-theta)
  logPrior <- log(theta) + 2*log(1-theta)
  return(logLikelihood + logPrior)
}
initVal <- 0.3
OptimRes <- optim(initVal,
                  logPostFn,
                  lower=0.2,upper = 0.9,
                  method=c("L-BFGS-B"),
                  control=list(fnscale=-1),#Maximize log posterior
                  hessian=TRUE)#Returns the hessian
postMode <- OptimRes$par
postSamples <- rnorm(n = 10000, mean = postMode,
                     sd = sqrt(solve(-OptimRes$hessian)))
lines(density(postSamples), col = 'red')


####2(b)####
#Plot the posterior distribution of θ 
#Draw from Gamma
shape <- sum(Wolves$y)+40
rate <- nrow(Wolves) + 2
thetas <- rgamma(10000, shape, rate)
hist(thetas)
#Compute the posterior probability that θ is smaller than 21
mean(thetas < 21)#0.1311


####2(c)####
#Do posterior inference in both models - to do posterior inference, we plot
#the postA vs postB
wolvesRegA <- Wolves[Wolves$x == 1,]
wolvesRegB <- Wolves[Wolves$x == 0,]

shape <- sum(wolvesRegA$y)+40
rate <- nrow(wolvesRegA) + 2
thetasRegA <- rgamma(10000, shape, rate)

shape <- sum(wolvesRegB$y)+40
rate <- nrow(wolvesRegB) + 2
thetasRegB <- rgamma(10000, shape, rate)

xlim <- c(min(thetasRegA), max(thetasRegB))
ylim <- c(0,1)
plot(density(thetasRegA))
plot(density(thetasRegA), xlim = xlim, ylim = ylim)
lines(density(thetasRegB), col = 'red')


# Calculate the probability that
# yB,j > yA,j on a randomly picked future week j.
#Draw posterior theta from region A and get y value for each theta 
predYRegA <- sapply(thetasRegA, function(x) rpois(1, x))
predYRegB <- sapply(thetasRegB, function(x) rpois(1, x))
mean(predYRegB > predYRegA)#0.6821


####2(d)####
mean(thetasRegB >= (thetasRegA*1.1))#0.9845


####3(a)####
#simulate 10000 draws from the joint posterior distribution
#of all regression coefficients and the error variance
mu_0 <- rep(0,14)
Omega_0 <- (1/(10^2))*diag(14)
v_0 <- 1
sigma2_0 <- 4^2#covariance
nIter <- 10000
postList <- BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter)
postVars <- postList$sigma2Sample
postBetas <- postList$betaSample
#Summarize the posterior of the regression coefficients by the point estimate 
#under the linear loss function - Median
apply(postBetas, 2, median)
#95% equal tail credible intervals for all parameters in β
apply(postBetas, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
#2nd column

####3(b)####
#posterior mean and posterior median of the standard deviation σ
mean(sqrt(postVars))
median(sqrt(postVars))

####3(c)####
#Create the new X
newCrim <- seq(min(X[,2]), max(X[,2]), 0.1)
n <- length(newCrim)
newX <- matrix(nrow = n, ncol = ncol(X))
for (i in 1:ncol(X)) {
  newX[,i] <- rep(XNewHouse[i], n)
}
newX[,2] <- newCrim

expValueMat <- postBetas %*% t(newX)
expValue95CredInt <- apply(expValueMat, 2, 
                           function(x) quantile(x, probs = c(0.025, 0.975)))
plot(x = newCrim, y = expValue95CredInt[1,], type = 'l',
     ylim = c(min(expValue95CredInt[1,]), max(expValue95CredInt[2,])),
     xlab = "Crime Rate",
     ylab = "Expected value of mu")
lines(x = newCrim, y = expValue95CredInt[2,])

####3(d)####
expValueMat <- as.vector(postBetas %*% as.matrix(XNewHouse, ncol = 1))
newSellingPrice <- sapply(1:nIter, 
                          function(x) rnorm(1, 
                                            expValueMat[i], 
                                            sqrt(postVars[i])))
mean(newSellingPrice >= 20)#0.9581

####3(e)####
#Use the posterior(beta,var) and x to get 890 y for each 10k combination
expValueMat <- postBetas %*% t(X)
yMat <- matrix(0, nrow = nIter, ncol = length(y))
for (i in 1:length(y)) {
  yMat[,i] <- sapply(1:nIter, 
                            function(x) rnorm(1,
                                              expValueMat[x, i],
                                              sqrt(postVars[x])))
  
}
tRep <- apply(yMat, 1, max)
mean(tRep > tAct)#0.4693
