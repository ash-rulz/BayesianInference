####1(b)####

post_pred_sim <- function(){
  #Step 1: Get theta from Gamma(2326, 7)
  theta <- rgamma(1, 2326,7)
  #Step 2: Get q from Pois(theta)
  q <- rpois(1, theta)
  return(q)
}
n <- 10000
q <- sapply(1:n, function(x) post_pred_sim())
#Step 3: Plot histogram
hist(q)

#Step 4: Pr(Q6 > 350|q1, ..., q5)
pr_q_gt350 <- sum(q > 350)/n
pr_q_gt350#0.1762

####1(c)####
a_seq <- seq(300, 500, 1)
utility_seq <- sapply(a_seq, function(x) mean(utility_func(x, q)))
plot(utility_seq, type = 'l', x = a_seq)
a_seq[which.max(utility_seq)]#363
abline(v = a_seq[which.max(utility_seq)])




####2(a)####
#Simulate draws from posterior
mu_0 <- as.vector(rep(0,6))
Omega_0 <- (1/100)*diag(6)
v_0 <- 1
sigma2_0 <- 100^2
nIter <- 10000
#Joint posterior draws for all betas and sigma2
joint_post_sim <- BayesLinReg(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter)

#Posterior mean for all parameters in beta
mean_betas <- apply(joint_post_sim$betaSample, 2, mean)
#-488.22719209   10.72500972   -0.02607327   33.22065501   -0.60755907   -0.13438897

#99% equal tail credible intervals for all parameters in beta
apply(joint_post_sim$betaSample, 2, function(x) quantile(x, 
                                                         probs = c(0.005, 
                                                                   0.995)))
####2(b) - Posterior mean and median of sigma2####
post_sd <- sqrt(joint_post_sim$sigma2Sample)
mean(post_sd)#39.96048
median(post_sd)#39.47453


####2(c)####
X1 <- seq(min(X[,2]), max(X[,2]), by = 0.1)
X0 <- rep(1, length(X1))
X1_sq <- X1^2
X2 <- rep(27, length(X1))
X2_sq <- X2^2
X1.X2 <- X1*X2
newX <- rbind(X0, X1, X1_sq,
              X2, X2_sq,
              X1.X2)
newMu_mat <- joint_post_sim$betaSample %*% newX
interval_mat <- apply(newMu_mat, 2, function(x) quantile(x, c(0.025, 0.975)))
plot(x = X1, y = interval_mat[2,], type = 'l')
lines(x = X1, y = interval_mat[1,], type = 'l')


####2(e)####
#Sim posterior predictive distribution for y
newX <- matrix(data = c(1, 50, 50^2,
               25, 25^2, 50*25), ncol = 1)
newY <- (joint_post_sim$betaSample %*% newX) +
  sapply(post_sd, function(x) rnorm(1, 0, x))
plot(newY)
hist(newY)

####2(f)####
#Calculate posterior predictive p-value
tAct <- max(y)
newMu <- joint_post_sim$betaSample %*% t(X)
tSim <- c()
for (i in 1:nIter) {
 tSim <- c(tSim, 
           max(rnorm(length(y), newMu[i, ], post_sd[i]))) 
}
mean(tSim >= tAct)#0.9901






####3(c)####
#Write a function in R that computes the log posterior distribution of θ
theta_post_dist1 <- function(theta){
  x <- c(.8, 1.1, .8, .9, 1)
  n <- length(x)
  logLikelihood <- n*log(theta) - theta * sum(x^3)
  logPrior <- 2 * log(theta) - 4 * theta
  return(exp(logLikelihood + logPrior))
}
thetas <- seq(0.1, 3, by = 0.01)
postDensity <- sapply(thetas, function(theta) theta_post_dist1(theta))

#Normalize the posterior density
normConst <- integrate(theta_post_dist1, lower = min(thetas), 
                        upper = max(thetas))
normPostDensity <- postDensity/normConst$value
plot(thetas, normPostDensity, type = 'l', col = 'blue')

#Confirming if the area under the normalized curve is 1
f <- approxfun(thetas, normPostDensity)
tot_area <- integrate(f, min(thetas), max(thetas))
print(paste("CDF:", tot_area$value))


####3(d)####
LogPostFn <- function(theta){
  x <- c(.8, 1.1, .8, .9, 1)
  n <- length(x)
  logLikelihood <- n*log(theta) - theta * sum(x^3)
  logPrior <- 2 * log(theta) - 4 * theta
  return(logLikelihood + logPrior)
}
initVal <- 0.1
OptimRes <- optim(initVal,
                  LogPostFn,
                  lower=0.1,
                  gr=NULL,
                  method=c("L-BFGS-B"),
                  control=list(fnscale=-1),#Maximize log posterior
                  hessian=TRUE)#Returns the hessian
postMode <- OptimRes$par
postVar <- solve(-OptimRes$hessian)
post_thetas <- rnorm(n = 4000,
                     postMode,
                     sqrt(solve(-OptimRes$hessian))
                     )
lines(density(post_thetas), col = 'red')
