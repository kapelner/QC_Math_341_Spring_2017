library(MCMCpack)
X = read.csv("sp_tot_ret_price_1950.csv")
head(X)
n = nrow(X)
x = na.omit(X[,4]) *100
graphics.off()

par(mfrow = c(4, 2))
hist(x, br = 10000, ylim = c(0, 30), xlim = c(-22, 13),
     main = "True data", xlab = "daily return (%)")
s = sd(x)
xbar = mean(x)

for (r in 1 : 7){
  theta = rnorm(1, xbar, s / sqrt(n))
  sigsq = rinvgamma(1, (n - 1) / 2, (n - 1) * s^2 / 2)
  xrep = rnorm(n, theta, sigsq)
  hist(xrep, br = 10000, ylim = c(0, 30), xlim = c(-22, 13), 
       main = paste("Replicate", r), xlab = "daily return (%)")
  
}

