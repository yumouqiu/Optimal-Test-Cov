library(MASS)
source("functions.R")

n = 150
p = c(1500)
signal.strength = seq(0.15, 0.30, 0.03)
rep = 1000
m = 0

QC.power = data.frame()
CJ.power = data.frame()
MTT.power = data.frame()

for (i1 in 1 : length(p)){
  p0 = p[i1]
  sigma02 = runif(p0, min = 1, max = 5)
  signal.sparsity = c(p0 / 2, p0 / 5, p0 / 10)
  for(i in 1 : length(signal.strength)){
    for (j in 1 : length(signal.sparsity)){
      Sigma = Cov.alternative.rho(p0, sigma02, signal.strength[i], signal.sparsity[j])
      eigen.res = eigen(Sigma)
      Sigma.half = eigen.res$vectors %*% diag(sqrt(eigen.res$values)) %*% solve(eigen.res$vectors)
      QC.res = c()
      CJ.res = c()
      MTT.res = c()
      for (k in 1 : rep){
        X0 = runif(n * p0, min = -1, max = 1)
        X0 = X0 * sqrt(3)
        X = matrix(X0, n, p0) %*% Sigma.half
        QC.res[k] = QC(X)
        CJ.res[k] = CJ(X)
        MTT.res[k] = MTT(X, 0.5)
      }
      m = m + 1
      QC.power[m, 1 : 4] = c(mean(QC.res), p0, signal.strength[i], j)
      CJ.power[m, 1 : 4] = c(mean(CJ.res), p0, signal.strength[i], j)
      MTT.power[m, 1 : 4] = c(mean(MTT.res), p0, signal.strength[i], j)
      
      cat("p = ", p0, "strength = ", signal.strength[i], "sparsity = ", j, "\n")
    }
  }
}

names(QC.power) = c("Power", "p", "Strength", "Sparsity")
names(CJ.power) = c("Power", "p", "Strength", "Sparsity")
names(MTT.power) = c("Power", "p", "Strength", "Sparsity")
QC.power$method = "QC"
CJ.power$method = "CJ"
MTT.power$method = "MTT"

Res = rbind(QC.power, CJ.power, MTT.power)
write.csv(Res, "power-uniform-n150-p1500.csv", row.names = FALSE)
