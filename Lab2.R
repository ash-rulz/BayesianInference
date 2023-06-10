####Question1: Linear Regression ####
library(mvtnorm)
library(readxl)
#Given - starting values for prior parameters
prior_precision <- 4.5 * diag(3)#omega_Not
prior_df <- 3#new_not
prior_sd0_sq <- 12#sd_not
#prior_mean <- c(-13, 115, -110)
prior_mean <- t(matrix(data = c(-7.532457, 
                                83.578936, 
                                -78.298468), nrow = 1))
#Given - data 
data <- readxl::read_excel("Linkoping2022.xlsx")
n <- nrow(data)
k <- 3+1 #3covariates + 1 intercept
data$time <- as.numeric(format(as.Date(data$datetime), "%j"))/365

#Part a - Get appropriate prior values that matches with the original regression
#Step 1: Get beta_0, beta_1, beta_2 by drawing from the posterior of beta(L3, S5)
#Step 1.a: draws from chi squared distribution
n_sim <- 50

plot_df <- as.data.frame(matrix(data = 0, nrow = 0, ncol = 4))
colnames(plot_df) <- c('prior_b0', 'prior_b1', 'prior_b2', 'y_pred')

for (i in 1:n_sim) {
  X <- rchisq(1, df = prior_df)
  #Step 1.b: compute prior_sd_sq
  prior_sd_sq <- (prior_df * prior_sd0_sq)/X
  #Step 1.c: plug prior_sd_sq into rmvnorm to get draws for betas
  prior_beta <- mvtnorm::rmvnorm(n = 1, mean = prior_mean,
                                 sigma = prior_sd_sq * solve(prior_precision))
  #print(prior_beta)
  #Step2: Get temp values from prior betas
  temp <- list(prior_beta[1] + prior_beta[2] * data$time + prior_beta[3] * data$time^2)
  plot_df <- rbind(plot_df, data.frame(
    prior_b0 = prior_beta[1], 
    prior_b1 = prior_beta[2], 
    prior_b2 = prior_beta[3],
    y_pred = I(temp)
  )
  )
}
#Plot the response data for prior betas 
plot(data$time, data$temp, pch = 20, xlab = "Time", ylab = "Temperature")
colors <- rainbow(nrow(plot_df))
for (i in 1:n_sim) {
  lines(data$time, plot_df$y_pred[[i]], col = colors[i], lwd = 0.5)
}
# legend("topright", legend = paste0("b0=", round(plot_df$prior_b0, 2),
#                                    ", b1=", round(plot_df$prior_b1, 2),
#                                    ", b2=", round(plot_df$prior_b2, 2)),
#        col = colors, lwd = 0.5, cex = 0.3)

#Part b(i) - Function to simulate draws from joint posterior distribution of 
#beta0, beta1, beta2, sd_sq & plot the marginal posterior of the parameters

#Task1: Function to simulate draws from joint posterior distribution of 
#beta0, beta1, beta2, sd_sq 
#Step 1: See L5, S5 to calculate the posterior parameters
#Calculate X/covariate matrix
covariate <- as.matrix(cbind(
  rep(1,nrow(data)),#Column 1 with intercept
  data$time,#time
  data$time^2#time squared
))
dim(covariate)
#Calculate beta_hat/OLSE
beta_hat <- solve(t(covariate)%*%covariate) %*% (t(covariate)%*%data$temp) 
dim(beta_hat)
#This is Linear Regression with conjugate prior
#Calculate posterior mean, precision, df, si
post_mean <- solve((t(covariate)%*%covariate) + prior_precision) %*%
  ((t(covariate)%*%covariate%*%beta_hat) + prior_precision%*%prior_mean)
post_precision <- (t(covariate) %*% covariate) + prior_precision
post_df <- prior_df+n
post_sdn_sq <- (prior_df * prior_sd0_sq +
                  ((t(data$temp)%*%data$temp) + 
                     (t(prior_mean)%*%prior_precision%*%prior_mean)-
                     t(post_mean)%*%post_precision%*%post_mean))/post_df
post_mean
post_sdn_sq
#The mean and the sd would be similar to the prior and data values
draw_from_join_post_fn <- function(){
  #Step 2: Draw from beta0, beta1, beta2, sd_sq
  #Step 2.a: draws from chi squared distribution
  X <- rchisq(1, df = post_df)
  #Step 2.b: compute post_sd_sq
  #Step 2.b: compute prior_sd_sq
  post_sd_sq <- as.numeric((post_df * post_sdn_sq)/X)
  #Step 2.c: plug post_sd_sq into rmvnorm to get draws for betas
  post_beta <- mvtnorm::rmvnorm(n = 1, mean = post_mean,
                                sigma = (post_sd_sq) * solve(post_precision))
  return(cbind(post_beta, post_sd_sq))
}
joint_post_df <- as.data.frame(matrix(nrow = 0, ncol = 4))
joint_post_df <- rbind(joint_post_df,
                       t(sapply(1:1000, function(x) draw_from_join_post_fn())))
joint_post_df#We get the posterior joint parameters beta0, beta1, beta2, sd_sq
colnames(joint_post_df) <- c('beta_0', 'beta_1', 'beta_2', 'sd')

#Task2: plot the marginal posterior of the parameters - L5, S4 - t-distribution
#This is for non-informative prior. For informative prior, 
#use the joint posterior itself.
s_sq <-  (t(data$temp - (covariate %*%beta_hat)) %*% 
  (data$temp - (covariate %*%beta_hat)))/(n-k)
s_sq <- as.numeric(s_sq)
#Marginal posterior beta drawn from the t-distribution
# mar_pos_beta <- mvtnorm::rmvt(10000,
#                               sigma = s_sq * solve(t(covariate)%*%covariate),
#                               df = n-k,
#                               delta = beta_hat)
hist(joint_post_df$beta_0, xlab = 'beta-0',
     main = 'beta-0 marginal posterior histogram')
hist(joint_post_df$beta_1, xlab = 'beta-1', 
     main = 'beta-1 marginal posterior histogram')
hist(joint_post_df$beta_2, xlab = 'beta-2', 
     main = 'beta-2 marginal posterior histogram')
hist(as.numeric(joint_post_df$sd), xlab = 'sd', 
     main = 'sd marginal posterior histogram', prob = TRUE)#Same as joint posterior distr
lines(density(joint_post_df$sd), col = 'blue')

#Part b(ii) - Plot the posterior median of the regression with the 5% and 95%
#posterior percentiles
#Marginal posterior beta %*% covariates gives the posterior response
pos_responses <- as.matrix(joint_post_df[,1:3]) %*% t(covariate)#1column for each day/x
pos_response_median <- apply(pos_responses, 2, median)
plot(data$time, data$temp, pch = 20, xlab = "Time", ylab = "Temperature")
lines(data$time, pos_response_median, col = 'red')#Median line
#Get the 90% equi-tailed posterior probability
pos_responses_ci = apply(pos_responses, 2, quantile, probs = c(0.05, 0.95))
lines(data$time, pos_responses_ci[1,], col = 'blue', lty = 2,)#5% CI
lines(data$time, pos_responses_ci[2,], col = 'blue', lty = 2,)#95% CI

#Part c: Find the time(covariate) which gives the maximum temp(response)
calc_x <- function(beta) {
  #Get the formula by taking the derivative of the reg expression w.r.t covariate
  #and equating it to 0 - optimization problem
  return(-beta[2] / (2 * beta[3]))
}
xTilde <- apply(joint_post_df[,1:3], 1, calc_x)
hist(xTilde, breaks = 25)


####Classification with logistic regression####
#Given
data <- read.table("WomenAtWork.dat", header = TRUE)
n <- nrow(data)
X <- as.matrix(data[,2:ncol(data)])
y <- data[,1]
dim(X)
Xnames <- colnames(X)
Npar <- dim(X)[2]
# Setting up the prior
mu <- as.matrix(rep(0,Npar)) # Prior mean vector
tau <- 2
Sigma <- (tau^2)*diag(Npar) # Prior covariance matrix

#Part a(i): Draw samples from posterior distribution - L6,S5
#Step1: Calculate the mean, beta_hat using numerical optimization - optim
# Functions that returns the log posterior for the logistic and probit regression.
# First input argument of this function must be the parameters we optimize on, 
# i.e. the regression coefficients beta.
LogPostLogistic <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas
  logLik <- sum( linPred*y - log(1 + exp(linPred)) )
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
}
# Select the initial values for beta
initVal <- matrix(0,Npar,1)
# The argument control is a list of options to the optimizer optim, 
# where fnscale=-1 means that we minimize 
# the negative log posterior. Hence, we maximize the log posterior.  
OptimRes <- optim(initVal,
                  LogPostLogistic,
                  gr=NULL,y, X, mu, Sigma,
                  method=c("BFGS"),
                  control=list(fnscale=-1),#Maximize log posterior
                  hessian=TRUE)#Returns the hessian
beta_hat <- OptimRes$par
print('The posterior mode is:')
print(beta_hat)
approxPostStd <- sqrt(diag(solve(-OptimRes$hessian))) # Computing approximate standard deviations.
print('The approximate posterior standard deviation is:')
print(approxPostStd)
print('Values of inverse of observed information at mode')
solve(-OptimRes$hessian)

post_samples <- mvtnorm::rmvnorm(n = 1000, 
                                 mean = beta_hat, 
                                 sigma = solve(-OptimRes$hessian)
                                 )
dim(post_samples)
# Part a(ii):Compute an approximate 95% equal tail posterior probability interval
# for the regression coefficient to the variable NSmallChild
post_NSmallChild_CI <- quantile(post_samples[,6], probs = c(0.05, 0.95))
post_NSmallChild_CI
glmModel <- glm(Work ~ 0 + ., data = data, family = binomial)
abs(glmModel$coefficients - apply(post_samples, 2, mean))

# Part b: Write a function that simulate draws from the  posterior predictive 
# distribution of Pr(y = 0|x), where the values of x are given
prob_not_working <- function(post_samples){
  new_X <- matrix(c(1, 18, 11, 7, 40, 1, 1), ncol = 1)
  linPred <- post_samples%*%new_X
  post_prob_working <- exp(linPred)/(1+exp(linPred))
  post_prob_not_working <- (1-post_prob_working)
  return(post_prob_not_working)
}
post_prob_not_working <- prob_not_working(post_samples)
hist(post_prob_not_working, main='Histogram of Pr(y = 0|x)')

# Part c: plot the posterior predictive distribution for the number of women, 
# out of these 13, that are not working
post_no_not_working_fn <- function(post_prob_not_working){
  post_no_not_working <- rbinom(n = length(post_prob_not_working), 
                                size = 13, prob = post_prob_not_working) 
  return(post_no_not_working)
}
post_13_not_working <- post_no_not_working_fn(post_prob_not_working)
hist(post_13_not_working)
