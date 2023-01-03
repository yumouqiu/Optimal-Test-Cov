library(MASS)
source("functions.R")

n = 150
p = c(500, 1000)
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
      QC.res = c()
      CJ.res = c()
      MTT.res = c()
      for (k in 1 : rep){
        X = mvrnorm(n, rep(0, p0), Sigma)
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
write.csv(Res, "power-normal-n150.csv", row.names = FALSE)
