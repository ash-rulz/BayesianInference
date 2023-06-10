####1(b)####
#Simulate 10000 draws from the predictive distribution of the maximal weight
nIter <- 10000
#Draw theta for 6th week
postTheta <- rgamma(1, 11, 11.25)
#Sim 10k draws of y for 6th week
predY6 <- rgamma(nIter, 2, postTheta)
mean(predY6 < 1.9)#0.5392

####1(c)####
nWeeks <- 30
postY30Mat <- matrix(0, nrow = nIter, ncol = nWeeks)
for (i in 1:nIter) {
  #Get posterior for future 30 weeks
  postThetas <- rgamma(nWeeks, 11, 11.25)
  #Get posterior prediction for each week
  postY30Weeks <- sapply(postThetas, function(theta) rgamma(1, 2, theta))  
  postY30Mat[i,] <- postY30Weeks
}
mean(apply(postY30Mat > 2.4, 1, sum))#10.5211


####1(d)####
cost <- seq(1,10,0.1)#Cost of building the escalator
loss <- c()
for (i in 1:length(cost)) {
  maxWeight <- 0.9 * log(cost[i])#Weight that can be held by the escalator
  #Number of weeks the escalator breaks out of the 30 weeks
  nBreakWeeks <- mean(apply(postY30Mat, 1, function(x) sum(x>maxWeight)))
  #Loss for the given cost
  loss <- c(loss, cost[i] + nBreakWeeks)
}
plot(cost, loss, type = 'l')
cost[which.min(loss)]#6


####2(a)####
#Simulate 10000 draws from the joint posterior
mu_0 <- as.vector(rep(0,8))
Omega_0 <- (1/9)*diag(8)
v_0 <- 1
sigma2_0 <- 9 
nIter <- 10000
X <- as.matrix(X)
postDraws <- BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter)
postBetas <- postDraws$betaSample
postVars <- postDraws$sigma2Sample
#Compute the 99% equal tail credible interval for β1 and interpret it
quantile(postBetas[,2], probs = c(0.5/100, 99.5/100))

####2(b)####
x1 <- x2 <- x5 <- const <-   1
x4 <- 0#Apartment sold in the south side, not inner city
x3 <- 0.5
newX <- matrix(c(const, x1, x2, x3, x4, x5, x1*x4, x1 * x5), ncol = 1)
newMu <- as.vector(postBetas %*% newX)
cv <- sqrt(postVars) /newMu 
median(cv)#0.4819864

####2(c)####
#Mu difference for inner city vs south side
const <- 1
x1 <- 1
x2 <- x3 <- 0
#Inner city
x4 <- 1
x5 <- 0#Not in the south side
newX <- matrix(c(const, x1, x2, x3, x4, x5, x1*x4, x1 * x5), ncol = 1)
newMuInnerCity <- as.vector(postBetas %*% newX)
#South side
x4 <- 0#Not in the inner city
x5 <- 1#In the south side
newX <- matrix(c(const, x1, x2, x3, x4, x5, x1*x4, x1 * x5), ncol = 1)
newMuSouth <- as.vector(postBetas %*% newX)
plot(density(newMuInnerCity), col ='blue', type = 'l', 
     xlim = c(min(newMuSouth),max(newMuInnerCity)))
lines(density(newMuSouth), col ='red', type = 'l')

#Effect of y from x1 is diff for south side vs !(south or inner)
const <- 1
x1 <- 1
x2 <- x3 <- 0#effect of x1 is being investigated
#South side
x4 <- 0#Not in the inner city
x5 <- 1#In the south side
newX <- matrix(c(const, x1, x2, x3, x4, x5, x1*x4, x1 * x5), ncol = 1)
newYSouth <- as.vector(postBetas %*% newX) + 
  sapply(postVars, function(var) rnorm(1, 0, sqrt(var)))
#!(South or Inner city)
x4 <- 0#Not in the inner city
x5 <- 0#Not on the south side
newX <- matrix(c(const, x1, x2, x3, x4, x5, x1*x4, x1 * x5), ncol = 1)
newYNeither <- as.vector(postBetas %*% newX) + 
  sapply(postVars, function(var) rnorm(1, 0, sqrt(var)))
plot(density(newYSouth), col ='blue', type = 'l', 
     xlim = c(min(newYNeither),max(newYSouth)))
lines(density(newYNeither), col ='red', type = 'l')
mean(newYSouth)
mean(newYNeither)

####2(d)####
const <- 1
x1 <- x2 <- -0.5
x3 <- 0
x4 <- 0#Not in the inner city
x5 <- 1#In the south side
newX <- matrix(c(const, x1, x2, x3, x4, x5, x1*x4, x1 * x5), ncol = 1)
newMuSouth <- as.vector(postBetas %*% newX)
hist(newMuSouth)
mean(newMuSouth > 0)


####2(e)####
x1 <- seq(min(X[,2]), max(X[,2]), by = 0.01)
n <- length(x1)
const <- rep(1, n)
x2 <- rep(1, n)
x3 <- rep(0.5, n)
x4 <- rep(1, n)#Inner city
x5 <- rep(0, n)#Not in the south side
newX <- as.matrix(cbind(const, x1, x2, x3, x4, x5, x1*x4, x1 * x5))
newYInner <- postBetas %*% t(newX)
for (i in 1:ncol(newYInner)) {
  newYInner[,i] <- newYInner[,i] + 
    sapply(postVars, function(var) rnorm(1, 0, sqrt(var)))
}
newYInnerInterval <- apply(newYInner, 2, 
                           function(x) quantile(x, probs = c(2.5/100, 97.5/100)))
plot(x1, newYInnerInterval[1,], type = 'l', col = 'blue', 
     ylim = c(min(newYInnerInterval[1,]), max(newYInnerInterval[2,])))
lines(x1, newYInnerInterval[2,], type = 'l', col = 'red')
