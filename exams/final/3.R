n = 24 * 60

true_lambda_calm = 5
true_lambda_busy = 5.6
true_m_calm_to_busy = 8 * 60
true_m_busy_to_calm = 19 * 60

x = c(
  rpois(true_m_calm_to_busy, true_lambda_calm),
  rpois(true_m_busy_to_calm - true_m_calm_to_busy, true_lambda_busy),
  rpois(n - true_m_busy_to_calm, true_lambda_calm)
)
graphics.off()
plot((1 : n), x, xlab = "time (in hours from midnight)", ylab = "calls per minute")


## hyperparams
alpha = 0
beta = 0

#chains
S = 1000
lambdacs = array(NA, S)
lambdabs = array(NA, S)
m1s = array(NA, S)
m2s = array(NA, S)
#start positions
lambdacs[1] = 3
lambdabs[1] = 4
m1s[1] = 10
m2s[1] = 1000

for (t in 2 : S){
  cat ("t = ", t, "\n")
  m1 = m1s[t - 1]
  m2 = m2s[t - 1]
  
  a = sum(x[1 : m1])
  b = sum(x[(m1 + 1) : m2])
  c = sum(x[(m2 + 1) : n])
  lambdac = rgamma(1, alpha + a + c, m1 + n - m2 + beta)
  lambdab = rgamma(1, alpha + b, m2 - m1)
  lambdacs[t] = lambdac
  lambdabs[t] = lambdab
  if (is.na(lambdac)){
    stop("lambdacbusted")
  }
  if (is.na(lambdab)){
    stop("lambdab busted")
  }
  
  #now we need to calculate all the m dist
  m1_dist = array(NA, n)
  const = sum(lgamma(x[1 : n] + 1))
  for (i in 1 : (m2 - 1)){
    a = sum(x[1 : i])
    b = sum(x[(i + 1) : m2])
    m1_dist[i] = a * log(lambdac) + b * log(lambdab) - lambdac * i  - lambdab * (m2 - i)
  }
  m1_dist = m1_dist - max(m1_dist, na.rm = TRUE)
  m1_dist = exp(m1_dist)
  m1_dist = m1_dist / sum(m1_dist, na.rm = TRUE)
  m1_dist[is.na(m1_dist)] = 0
  # plot(m1_dist)
  
  m1 = sample(1 : n, 1, prob = m1_dist)
  m1s[t] = m1
  
  m2_dist = array(NA, n)
  for (i in (m1 + 1) : (n - 1)){
    b = sum(x[(m1 + 1) : i])
    c = sum(x[(i + 1) : n])
    m2_dist[i] = b * log(lambdab) +  c * log(lambdac) - lambdab * i - lambdac * (n - i)
  }
  m2_dist = m2_dist - max(m2_dist, na.rm = TRUE)
  m2_dist = exp(m2_dist)
  m2_dist = m2_dist / sum(m2_dist, na.rm = TRUE)
  m2_dist[is.na(m2_dist)] = 0
  # plot(m2_dist)
  
  m2s[t] = sample(1 : n, 1, prob = m2_dist)
  
  if (m2s[t] < m1){
    stop("m2 busted")
  }
}



###assess convergence
par(mfrow = c(4, 1))
S0 = 1
Sf = 80
plot(S0 : Sf, lambdacs[S0 : Sf], xlab = "iteration #", ylab = "lambda_c")
# abline(h = mean(lambdacs[B : S0]), col = "blue")
# abline(h = true_lambda_1, col = "red")
# abline(v = B, col = "grey")

plot(S0 : Sf, lambdabs[S0 : Sf], xlab = "iteration #", ylab = "lambda_b")
# abline(h = mean(lambdabs[B : S0]), col = "blue")
# abline(h = true_lambda_2, col = "red")
# abline(v = B, col = "grey")

plot(S0 : Sf, m1s[S0 : Sf], xlab = "iteration #", ylab = "m_1")
# abline(h = mean(ms[B : S0]), col = "blue")
# abline(h = sqrt(true_m), col = "red")
# abline(v = B, col = "grey")

plot(S0 : Sf, m2s[S0 : Sf], xlab = "iteration #", ylab = "m_2")
# abline(h = mean(ms[B : S0]), col = "blue")
# abline(h = sqrt(true_m), col = "red")
# abline(v = B, col = "grey")
#plot
B = 10

##assess autocorrelation

par(mfrow = c(4, 1))
K = 5
acf(lambdacs[B : S], lag.max = K, xlim = c(1, K))
acf(lambdabs[B : S], lag.max = K, xlim = c(1, K))
acf(m1s[B : S], lag.max = K, xlim = c(1, K))
acf(m2s[B : S], lag.max = K, xlim = c(1, K))
T = 3
#plot

#burn and thin
lambdacs = lambdacs[B : S]
lambdacs = lambdacs[seq(1, S - B, by = T)]
lambdabs = lambdabs[B : S]
lambdabs = lambdabs[seq(1, S - B, by = T)]
m1s = m1s[B : S]
m1s = m1s[seq(1, S - B, by = T)]
m2s = m2s[B : S]
m2s = m2s[seq(1, S - B, by = T)]


#look at posteriors with post-exp at 95% CI
par(mfrow = c(4, 1))

breaks = 300
hist(lambdacs, br = breaks)
abline(v = mean(lambdacs), col = "blue", lwd = 3)
abline(v = quantile(lambdacs, 0.025), col = "grey", lwd = 3)
abline(v = quantile(lambdacs, 0.975), col = "grey", lwd = 3)
# abline(v = true_lambda_1, col = "red", lwd = 3)

hist(lambdabs, br = breaks)
abline(v = mean(lambdabs), col = "blue", lwd = 3)
abline(v = quantile(lambdabs, 0.025), col = "grey", lwd = 3)
abline(v = quantile(lambdabs, 0.975), col = "grey", lwd = 3)
# abline(v = true_lambda_2, col = "red", lwd = 3)

hist(m1s, br = breaks)
abline(v = mean(m1s), col = "blue", lwd = 3)
abline(v = quantile(m1s, 0.025), col = "grey", lwd = 3)
abline(v = quantile(m1s, 0.975), col = "grey", lwd = 3)
# abline(v = true_m, col = "red", lwd = 3)


hist(m2s, br = breaks)
abline(v = mean(m2s), col = "blue", lwd = 3)
abline(v = quantile(m2s, 0.025), col = "grey", lwd = 3)
abline(v = quantile(m2s, 0.975), col = "grey", lwd = 3)
#plot


