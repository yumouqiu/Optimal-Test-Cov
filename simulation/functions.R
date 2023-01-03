rho2.density = function(z, n){
  z2.density = z^2 * (1 - z^2)^((n - 4) / 2) 
  return(z2.density)
}

Threshold = function(X, s){
  # Thresholding statistic for normal data
  n = dim(X)[1]
  p = dim(X)[2]
  q = p * (p - 1) / 2
  lambda = 4 * s * log(p)
  sample.cor = cor(X)
  diag(sample.cor) = 0
  M = n * sample.cor^2
  Test.stat = sum(M * (M > lambda)) / 2
  lower.limit = sqrt(lambda / n)
  mu0 = 2 * n * q * integrate(rho2.density, lower = lower.limit, upper = 1, n = n)$value / beta(1 / 2, (n - 2) / 2)
  sigma20 = q * ( 2 * (lambda^(3 / 2) + 3 * lambda^(1 / 2)) * dnorm(sqrt(lambda)) + 6 * (1 - pnorm(sqrt(lambda))) )
  return((Test.stat - mu0) / sqrt(sigma20))
}

MTT = function(X, s0){
  # MTT statistic for normal data
  p = dim(X)[2]
  s1 = seq(s0, 0.99, 0.005)
  Single.Threshold = c()
  for(i in 1 : length(s1)){
    s = s1[i]
    Single.Threshold[i] = Threshold(X, s)
  }
  Test.stat = max(Single.Threshold)
  a1 = sqrt(2 * log(log(p)))
  b1 = 2 * log(log(p)) + 0.5 * log(log(log(p))) - 0.5 * log(pi) + log(1 - s0 - 0.01)
  reject = 1 * (Test.stat > (-log(log((1 - 0.05)^(-1))) + b1) / a1)
  return(reject)
}

Threshold.nonnormal = function(X, s){
  # Thresholding statistic for nonnormal data using sample correlations 
  # under the case that all variables are independent under the null hypothesis
  n = dim(X)[1]
  p = dim(X)[2]
  q = p * (p - 1) / 2
  lambda = 4 * s * log(p)
  sample.cor = cor(X)
  diag(sample.cor) = 0
  M = n * sample.cor^2
  Test.stat = sum(M * (M > lambda)) / 2
  mu0 = q * (2 * lambda^(1 / 2) * dnorm(sqrt(lambda)) + 2 * (1 - pnorm(sqrt(lambda))) )
  sigma20 = q * ( 2 * (lambda^(3 / 2) + 3 * lambda^(1 / 2)) * dnorm(sqrt(lambda)) + 6 * (1 - pnorm(sqrt(lambda))) )
  return((Test.stat - mu0) / sqrt(sigma20))
}

MTT.nonnormal = function(X, s0){
  p = dim(X)[2]
  s1 = seq(s0, 0.99, 0.005)
  Single.Threshold = c()
  for(i in 1 : length(s1)){
    s = s1[i]
    Single.Threshold[i] = Threshold.nonnormal(X, s)
  }
  Test.stat = max(Single.Threshold)
  a1 = sqrt(2 * log(log(p)))
  b1 = 2 * log(log(p)) + 0.5 * log(log(log(p))) - 0.5 * log(pi) + log(1 - s0 - 0.01)
  reject = 1 * (Test.stat > (-log(log((1 - 0.05)^(-1))) + b1) / a1)
  return(reject)
}

CJ = function(X){
  n = dim(X)[1]
  p = dim(X)[2]
  sample.cor = cor(X)
  diag(sample.cor) = 0
  Test.stat = max(abs(sample.cor))
  reject.value = n^(-1) * (4 * log(p) - log(log(p)) - log(8 * pi) - 2 * log(log((1 - 0.05)^(-1))))
  return(1 * (Test.stat^2 > reject.value))
}

QC = function(X){
  dyn.load("bandcov.so")
  source("bandcov.R")
  Test.stat = bandtest.stat(X, 0)
  reject = 1 * (Test.stat > 2 * qnorm(0.95))
  return(reject)
}

Cov.alternative.rho = function(p, Sigma, rho, number){
  Cov = diag(1, p)
  q = p * (p - 1) / 2
  for (i in 1 : number){
    Cov[i, i + 1] = rho
    Cov[i + 1, i] = rho
  }
  return(diag(sqrt(Sigma)) %*% Cov %*% diag(sqrt(Sigma)))
}

Cov.alternative.rho.1 = function(p, Sigma, rho, number){
  Cov = diag(1, p)
  q = p * (p - 1) / 2
  for (i in 1 : (p - 1)){
    Cov[i, i + 1] = rho
    Cov[i + 1, i] = rho
  }
  for (i in 1 : (p - 2)){
    Cov[i, i + 2] = rho
    Cov[i + 2, i] = rho
  }
  for (i in 1 : (number - (2 * p - 3))){
    Cov[i, i + 3] = rho
    Cov[i + 3, i] = rho
  }
  return(diag(sqrt(Sigma)) %*% Cov %*% diag(sqrt(Sigma)))
}

# Cov.alternative = function(p, r, beta){
#   Cov = diag(1.5, p)
#   q = p * (p - 1) / 2
#   strength = sqrt(4 * r * log(p) / n)
#   number = floor(q^(1 - beta))
#   for (i in 1 : number){
#     Cov[i, i + 1] = strength
#     Cov[i + 1, i] = strength
#   }
#   return(Cov)
# }
# 
# Cov.alternative1 = function(p = 200, r, beta){
#   Cov = diag(1.5, p)
#   q = p * (p - 1) / 2
#   strength = sqrt(4 * r * log(p) / n)
#   number = floor(q^(1 - beta))
#   for (i in 1 : (p - 1)){
#     Cov[i, i + 1] = strength
#     Cov[i + 1, i] = strength
#   }
#   for (i in 1 : (number - (p - 1))){
#     Cov[i, i + 2] = strength
#     Cov[i + 2, i] = strength
#   }
#   return(Cov)
# }
# 
# Cov.alternative2 = function(p = 400, r, beta){
#   Cov = diag(1.5, p)
#   q = p * (p - 1) / 2
#   strength = sqrt(4 * r * log(p) / n)
#   number = floor(q^(1 - beta))
#   for (i in 1 : (p - 1)){
#     Cov[i, i + 1] = strength
#     Cov[i + 1, i] = strength
#   }
#   for (i in 1 : (p - 2)){
#     Cov[i, i + 2] = strength
#     Cov[i + 2, i] = strength
#   }
#   for (i in 1 : (number - (2 * p - 3))){
#     Cov[i, i + 3] = strength
#     Cov[i + 3, i] = strength
#   }
#   return(Cov)
# }