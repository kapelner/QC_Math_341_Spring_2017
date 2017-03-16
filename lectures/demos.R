

lambda = .1

Nsim = 100000

x = array(NA, Nsim)


for (i in 1 : Nsim){
  beta = rexp(1, lambda)
  x[i] = rbeta(1, 1, beta)
}

hist(x, br = 1000)


Nsim = 100000

x = array(NA, Nsim)


for (i in 1 : Nsim){
  if (runif(1) < 1/3){
    x[i] = rnorm(1)
  } else {
    x[i] = rexp(1, 1/3)
  }
}

hist(x, br = 1000)