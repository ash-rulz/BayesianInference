#Question 1
#Bernoulli Data
n <- 70 #Total trials
s <- 22 #Successes
f <- n-s #Failures

#Prior
pr_alpha <- pr_beta <- 8

#Posterior
pos_alpha <- pr_alpha + s
pos_beta <- pr_beta + f

#1(a)
pos_mean <- c()
pos_sd <- c()
n_draws <- seq(from = 10, to = 10000, by = 10)
for (i in n_draws) {
  #Draw random var for posterior
  pos_rv <- rbeta(n = i, pos_alpha, pos_beta)
  pos_mean <- c(pos_mean, mean(pos_rv))
  pos_sd <- c(pos_sd, sd(pos_rv))
}
plot(n_draws, pos_mean, type = 'l',
     xlab = 'No of draws', ylab = 'Posterio mean')
true_mean <- pos_alpha/(pos_alpha+pos_beta)
abline(h = true_mean, col = 'red')
#The pos_mean converges to true_mean as draws increases

plot(n_draws, pos_sd, type = 'l',
     xlab = 'No of draws', ylab = 'Posterio standard deviation')
true_sd <- sqrt((pos_alpha * pos_beta)/((pos_alpha + pos_beta)^2 * (pos_alpha + pos_beta + 1)))
abline(h = true_sd, col = 'red')
#The pos_sd converges to true_mean as draws increases

#1(b)
pos_rv <- rbeta(n = 10000, pos_alpha, pos_beta)
#Filter out all theta's greater than 0.3
filt_pos_rv <- pos_rv[pos_rv>0.3]

#Getting the probability of theta>0.3
p_gt0.3 <- length(filt_pos_rv)/length(pos_rv)

#To get the exact/true value, using the pbeta to get the CDF.
cdf_0.3 <- pbeta(0.3, pos_alpha, pos_beta, lower.tail = FALSE)
print(paste('The prob > 0.3 is:', round(p_gt0.3,2)))
print(paste('CDF theta > 0.3 is:', round(cdf_0.3,2)))

#1(c)
pos_rv <- rbeta(n = 10000, pos_alpha, pos_beta)

#Calculating the phi values from theta values
phi_rv <- pos_rv/(1-pos_rv)

#Plotting the distribution of phi
hist(phi_rv, breaks = 40, xlab = "odds")

phi_den <- density(phi_rv)
plot(phi_den, col = 'red', 
     main = "Kernel density estimation of phi")#Question: Is this what is expected?

#Question 2
#2(a)
#Normal model with unknown variance
n_draws <- 10000
data_mean <- 3.6
Y <- c(33, 24, 48, 32, 55, 74, 23, 17)
n <- length(Y)

#Step1: draws from chi squared distribution
X <- rchisq(n_draws, df = n)

#Step2: Compute sigma^2 - draw from inverse chi-squared
taosq <- sum((log(Y)-data_mean)^2)/n
pos_sigma_sq_rv <- (n * taosq)/X

#2(b)
#phi is CDF for normal standard normal with 0 mean and unit variance - hence pnorm
phi <- pnorm(sqrt(pos_sigma_sq_rv/2), mean = 0, sd = 1)
G_rv <- 2 * phi - 1
hist(G_rv, breaks = 40)

#2(c)
cred_int <- quantile(G_rv, probs = c(0.025, 0.975))
hist(G_rv, breaks = 100)
abline(v = c(cred_int[1], cred_int[2]), col = 'red')

#2(d)
#Computing HPDI

#Getting the kernel density
ker_den <- density(G_rv)
plot(ker_den)

#Calculating the total area under the curve
f <- approxfun(ker_den$x, ker_den$y)
tot_area <- integrate(f, min(ker_den$x), max(ker_den$x))
print(tot_area)

#Sort the densities and the x-values based on descending order of density values
sorted_index <- order(ker_den$y, decreasing = TRUE)
ker_den_y <- ker_den$y[sorted_index]
ker_den_x <- ker_den$x[sorted_index]

#Finding the indexes which gives the area <= 0.95
cum_area <- 0
indx <- 0
while (cum_area <= .95) {
  f <- approxfun(ker_den_x, ker_den_y)
  indx <- indx + 1
  #Getting the area of the sorted pairs of indexes
  area <- integrate(f, ker_den_x[indx], ker_den_x[indx+1])
  cum_area <- area$value/tot_area$value
  #print(paste("Index:", indx, " cum_area:", cum_area))
}
print(paste("Lower limit:", ker_den_x[indx]))
print(paste("Upper limit:", ker_den_x[indx+1]))
hdpi <- HDInterval::hdi(ker_den, credMass = 0.95)
plot(ker_den)
abline(v = c(ker_den_x[indx], ker_den_x[indx+1]))
abline(v = c(hdpi[1], hdpi[2]), col = 'red')
abline(v = c(cred_int[1], cred_int[2]), col = 'blue')

#Question 3
#3(a)
#Posterior
pos_pdf <- function(k){
  #Data
  Y <- c(-2.79, 2.33, 1.83, -2.44, 2.23, 2.33, 2.07, 2.02, 2.14, 2.5)
  data_mean <- 2.4
  n <- length(Y)
  bes_val <- besselI(k, nu = 0)
  return((1/(bes_val)^n) * exp(k * sum(cos(Y-data_mean)) - (k/2)))
}
k = seq(0, 10, 0.01)
#Calculating the normalizing constant for given unknown distribution
norm_const <- integrate(pos_pdf, lower = min(k), upper = max(k))

#Getting the normalized pdf
norm_pos_pdf <- pos_pdf(k)/norm_const$value
plot(k, norm_pos_pdf, type = 'l')

#Confirming that the area is 1  
f <- approxfun(k, norm_pos_pdf)
tot_area <- integrate(f, min(k), max(k))
print(paste("CDF:", tot_area$value))

#3(b)
#The mode is the k with the maximum frequency/probability
print(paste("The mode:", k[which.max(norm_pos_pdf)]))


