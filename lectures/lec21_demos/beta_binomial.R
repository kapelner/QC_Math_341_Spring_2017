library(VGAM)
n = 100
x = 29
alpha = 1
beta = 1
plot(1:n, dbetabinom.ab(1:n, n, x + alpha, n-x + beta))