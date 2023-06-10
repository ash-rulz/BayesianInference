####1(a)####
x <- c(13, 8, 11, 7)
sum(x)
####1(d)###
logPostFn <- function(theta, x, N){
  n <- length(x)
  logLik <- sum(x)*log(theta) + 
    (n*N-sum(x))*log(1-theta)
  logPrior <- log(theta) + 2*log(1-theta)
  return(logLik + logPrior)
}
N <- 20
thetaGrid <- seq(0.2, 1, 0.01)
thetaPostDenUnnorm <- exp(logPostFn(thetaGrid, x, N))

thetaPostDenNorm <- thetaPostDenUnnorm/(0.01*sum(thetaPostDenUnnorm))
#Confirming that the area is 1
f <- approxfun(thetaGrid, thetaPostDenNorm)
tot_area <- integrate(f, min(thetaGrid), max(thetaGrid))
print(paste("CDF:", tot_area$value))

plot(thetaGrid, thetaPostDenNorm, type='l')

####1(e)####
initVal <- 0.1
OptimRes <- optim(initVal,
                  logPostFn,
                  gr=NULL, x, N, lower=0.2, upper = 0.99,
                  method=c("L-BFGS-B"),
                  control=list(fnscale=-1),#Maximize log posterior
                  hessian=TRUE)#Returns the hessian
thetaPostApprox <-  dnorm(thetaGrid, 
                          mean = OptimRes$par, 
                          sd = sqrt(solve(-OptimRes$hessian)))
lines(thetaGrid, thetaPostApprox, col = 'red')




####2(a)####
sum(Wolves$y)+40
dim(Wolves)

####2(b)####
thetaPost <- rgamma(10000, 3383, 158)
hist(thetaPost)
mean(thetaPost<21)#0.1349

####2(c)####

thetaPostA <- rgamma(10000, 
                     sum(Wolves[Wolves$x == 1,]$y)+40,
                     sum(Wolves$x == 1) + 2
)
thetaPostB <- rgamma(10000, 
                     sum(Wolves[Wolves$x == 0,]$y)+40,
                     sum(Wolves$x == 0) + 2
)
plot(density(thetaPostA),xlim=c(16, 26), ylim = c(0,1))
lines(density(thetaPostB), col = 'red')

yA <- rpois(10000, thetaPostA)
yB <- rpois(10000, thetaPostB)
mean(yB>yA)#6905

####2(d)####
mean(yB) > mean(yA+0.1*yA)
mean(yA+0.1*yA)


####3(a)####
#Simulate 10000 draws from the joint posterior distribution
mu_0 <- rep(0,14)
Omega_0 <- (1/10^2)*diag(14)
v_0 <- 1
sigma2_0 <- 4^2
nIter <- 10000
postResult <- BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter)
postbetaMat <- postResult$betaSample
postVarVect <- postResult$sigma2Sample
#Point estimate of posterior regression coefficient - median
apply(postbetaMat, 2, median)
#95% equal tail credible intervals for all parameters
apply(postbetaMat, 2, 
      function(beta) quantile(beta, probs = c(0.025, 0.975)))
####3(b)####
mean(sqrt(postVarVect))
median(sqrt(postVarVect))
####3(c)####
newCrim <- seq(min(X[,2]), max(X[,2]), 0.1)
const <- rep(1, length(newCrim))
newX <- cbind(const, newCrim)
restCov <- matrix(c(40,1.5,0,0.5,6,30,5,3,300,17,390,4),nrow = 1)
restCovMat <- apply(restCov, 2, function(x) rep(x,length(newCrim)))
newX <- cbind(newX, restCovMat)
dim(newX)
newMu <- postbetaMat %*% t(newX)
newMuInt <- apply(newMu, 2, 
      function(mu) quantile(mu, probs=c(0.025, 0.975)))
plot(newCrim, newMuInt[1,], type = 'l', ylim = c(12, 35))
lines(newCrim, newMuInt[2,], type = 'l')
####3(d)####
newMu <- postbetaMat%*%XNewHouse
newY <- rnorm(nIter, newMu, sqrt(postVarVect))
mean(newY>20)#0.9432

####3(e)####
newMu <- postbetaMat%*%t(X)
newY <- matrix(NA, nrow = nIter, ncol = ncol(newMu))
for (i in 1:ncol(newMu)) {
  newY[,i] <- sapply(1:nIter, 
                     function(j) rnorm(1, 
                                       newMu[j,i],
                                       sqrt(postVarVect[j])))
}
tYRep <- apply(newY, 1, max)
tYAct <- max(y)
mean(tYRep >= tYAct)#0.4715
