library(MASS)
source("functions.R")

n = c(100, 150)
p = c(1500)
sigma02 = list()
rep = 1000
j = 0

QC.power = data.frame()
CJ.power = data.frame()
MTT.power = data.frame()

for (i1 in 1 : length(p)){
  p0 = p[i1]
  sigma02[[1]] = rep(1, p0)
  sigma02[[2]] = runif(p0, min = 1, max = 5)
  sigma02[[3]] = runif(p0, min = 1, max = 10)
  for (i2 in 1 : length(n)){
  	for (i3 in 1 : length(sigma02)){
      n0 = n[i2]
      Sigma = diag(sigma02[[i3]])
      QC.res = c()
      CJ.res = c()
      MTT.res = c()
      
      for (k in 1 : rep){
        X0 = runif(n0 * p0, min = -1, max = 1)
        X0 = X0 * sqrt(3)
        X = matrix(X0, n0, p0) %*% diag(sqrt(sigma02[[i3]]))
        QC.res[k] = QC(X)
        CJ.res[k] = CJ(X)
        MTT.res[k] = MTT.nonnormal(X, 0.5)
      }
      
      j = j + 1
      QC.power[j, 1] = mean(QC.res)
      QC.power[j, 2] = n0
      QC.power[j, 3] = p0
      QC.power[j, 4] = i3
      
      CJ.power[j, 1] = mean(CJ.res)
      CJ.power[j, 2] = n0
      CJ.power[j, 3] = p0
      CJ.power[j, 4] = i3
      
      MTT.power[j, 1] = mean(MTT.res)
      MTT.power[j, 2] = n0
      MTT.power[j, 3] = p0
      MTT.power[j, 4] = i3
      
      cat("n = ", n0, "p = ", p0, "sigma2 = ", i3, "\n")
    }
  }
}

names(QC.power) = c("Size", "n", "p", "sigma2")
names(CJ.power) = c("Size", "n", "p", "sigma2")
names(MTT.power) = c("Size", "n", "p", "sigma2")
QC.power$method = "QC"
CJ.power$method = "CJ"
MTT.power$method = "MTT"

Res = rbind(QC.power, CJ.power, MTT.power)
write.csv(Res, "size-uniform-p1500.csv", row.names = FALSE)

