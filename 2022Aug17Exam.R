####1(d)####
#Function for log posterior distribution of θ as the sum
# of the log likelihood and the log prior distribution of θ
logPostFn <- function(theta){
  x <- c(0.7, 1.1, 0.9, 1.5)
  n <- length(x)
  logLikelihood <- 3*n*log(theta) - theta*sum(1/x) 
  logPrior <- 3*log(theta) - 2*theta
  return(exp(logLikelihood + logPrior))
}

thetas <- seq(0.1, 5, 0.01)

#Calculating the normalizing constant for given unknown distribution
norm_const <- integrate(logPostFn, lower = min(thetas), upper = max(thetas))

#Getting the normalized pdf
norm_pos_pdf <- logPostFn(thetas)/norm_const$value

plot(thetas, norm_pos_pdf, type = 'l', xlab = 'Theta'
     ,ylab = 'Normalized Posterior Density')

####1(e)####
logPostFn <- function(theta){
  x <- c(0.7, 1.1, 0.9, 1.5)
  n <- length(x)
  logLikelihood <- 3*n*log(theta) - theta*sum(1/x) 
  logPrior <- 3*log(theta) - 2*theta
  return(logLikelihood + logPrior)
}

initVal <- 0.1
OptimRes <- optim(initVal,
                  logPostFn,
                  gr=NULL,
                  lower=0.1,
                  method=c("L-BFGS-B"),
                  control=list(fnscale=-1),#Maximize log posterior
                  hessian=TRUE)#Returns the hessian
postMode <- OptimRes$par
postSamples <- rnorm(n = 10000, mean = postMode, 
                     sd = sqrt(solve(-OptimRes$hessian)))
lines(density(postSamples), col = 'red')


####2(a)####
#Simulate 20000 draws from the joint posterior
#Prior
mu_0 <- c(0,0,0)
Sigma_0 <- 10^2 * diag(3)

betaSamples <-  as.data.frame((BayesLogitReg(y, X, mu_0, Sigma_0, 20000))$betaSample)
colnames(betaSamples) <- c('beta0', 'beta1', 'beta2')

#Compute a 95% equal tail credible interval for β1 and interpret it.
hist(betaSamples$beta1)
quantile(betaSamples$beta1, probs = c(0.025,0.975))
#0.01336806 0.18295450 

####2(b)####
#Compute the joint posterior probability that both β1 > 0 and β2 > 0
nrow(betaSamples[(betaSamples$beta1>0) &
              (betaSamples$beta2>0),])/nrow(betaSamples)#0.91545

####2(c)####
newX <- matrix(data = c(1, 5, 1), ncol = 3)
betaSamples <- as.matrix(betaSamples)
linPred <- betaSamples%*%t(newX)
#Post prob that bridge i needs to be repaired within the next 5 years
postProb1 <- exp(linPred)/(1+exp(linPred))
#Post prob that bridge i need not be repaired within the next 5 years
postProb0 <- as.vector((1-postProb1)/postProb1)
hist(postProb0)
newX <- as.data.frame(cbind(X, postDistr))
newX <- newX[newX$x2 == 1,]
plot(x = newX$x1, y = newX$postDistr, xlab = "Bridge age",
     ylab = "Post Distribution")

####2(d)####
x1 <- seq(min(X[,2]), max(X[,2]), by = 0.1)
const <- rep(1, length(x1))
x2 <- rep(0, length(x1))
newX <- as.matrix(data.frame(const = const, x1 = x1, x2 = x2))

linPred <- betaSamples%*%t(newX)
postProb <- exp(linPred)/(1+exp(linPred))
postProbInt <- apply(postProb, 2, function(x) quantile(x, c(0.025, 0.975)))
plot(x = x1, y = postProbInt[2,], type = 'l', col = 'blue')
y_min <- min(postProbInt)
y_max <- max(postProbInt)
ylim <- c(y_min, y_max)
plot(x = x1, y = postProbInt[1,], type = 'l', col = 'blue', ylim = ylim)
lines(x = x1, y = postProbInt[2,], col = 'red')

####2(e)####
newX <- matrix(data = c(1, 40, 1), ncol = 3)
linPred <- betaSamples%*%t(newX)
postProb <- as.vector(exp(linPred)/(1+exp(linPred)))
hist(postProb)
mean(postProb > 0.5)#0.1527


####3(c)####