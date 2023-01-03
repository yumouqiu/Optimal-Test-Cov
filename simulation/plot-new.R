library(tidyverse)
library(ggthemes)
library(ggpubr)

A1 = read.csv("power-normal-n100.csv")
A2 = read.csv("power-normal-n150.csv")
A3 = read.csv("power-normal-n100-p1500.csv")
A4 = read.csv("power-normal-n150-p1500.csv")
A1$n = "n = 100"
A2$n = "n = 150"
A3$n = "n = 100"
A4$n = "n = 150"
A = rbind(A1, A2, A3, A4)

#A$p[A$p == 500] = "p = 500"
#A$p[A$p == 1000] = "p = 1000"
A$p = as.factor(A$p)
A$Sparsity[A$Sparsity == 1] = "p / 2"
A$Sparsity[A$Sparsity == 2] = "p / 5"
A$Sparsity[A$Sparsity == 3] = "p / 10"

A$method = factor(A$method, levels = c("MTT", "QC", "CJ"))
A$p = factor(A$p, levels = c("500", "1000", "1500"))
A$Sparsity = factor(A$Sparsity, levels = c("p / 2", "p / 5", "p / 10"))

p1 = A %>% 
  ggplot(aes(y = Power, x = Strength, col = method)) + 
  geom_line(aes(linetype = p)) + geom_point(aes(shape = method)) + 
  facet_grid(n~Sparsity, 
             labeller = labeller(.rows = label_value, .cols = label_both)) + 
  geom_hline(yintercept = 0.05, colour = "black", linetype = "dashed") + 
  ggtitle("Normal distribution") + xlab("Correlation (signal strength)") + 
  scale_x_continuous(breaks = c(0.15, 0.18, 0.21, 0.24, 0.27, 0.3)) + 
  theme_few() + theme(legend.position = "bottom") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12))



A1 = read.csv("power-uniform-n100.csv")
A2 = read.csv("power-uniform-n150.csv")
A3 = read.csv("power-uniform-n100-p1500.csv")
A4 = read.csv("power-uniform-n150-p1500.csv")
A1$n = "n = 100"
A2$n = "n = 150"
A3$n = "n = 100"
A4$n = "n = 150"
A = rbind(A1, A2, A3, A4)

#A$p[A$p == 500] = "p = 500"
#A$p[A$p == 1000] = "p = 1000"
A$p = as.factor(A$p)
A$Sparsity[A$Sparsity == 1] = "p / 2"
A$Sparsity[A$Sparsity == 2] = "p / 5"
A$Sparsity[A$Sparsity == 3] = "p / 10"

A$method = factor(A$method, levels = c("MTT", "QC", "CJ"))
A$p = factor(A$p, levels = c("500", "1000", "1500"))
A$Sparsity = factor(A$Sparsity, levels = c("p / 2", "p / 5", "p / 10"))

p2 = A %>% 
  ggplot(aes(y = Power, x = Strength, col = method)) + 
  geom_line(aes(linetype = p)) + geom_point(aes(shape = method)) + 
  facet_grid(n~Sparsity, 
             labeller = labeller(.rows = label_value, .cols = label_both)) + 
  geom_hline(yintercept = 0.05, colour = "black", linetype = "dashed") + 
  ggtitle("Uniform distribution") + xlab("Correlation (signal strength)") + 
  scale_x_continuous(breaks = c(0.15, 0.18, 0.21, 0.24, 0.27, 0.3)) + 
  theme_few() + theme(legend.position = "bottom") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Figure 2, normal distribution case
p1

# Figure 2, uniform distribution case
p2


A1 = read.csv("power-normal-n100-dense.csv")
A2 = read.csv("power-normal-n150-dense.csv")
A3 = read.csv("power-uniform-n100-dense.csv")
A4 = read.csv("power-uniform-n150-dense.csv")
A1$setting = "Normal, n = 100"
A2$setting = "Normal, n = 150"
A3$setting = "Uniform, n = 100"
A4$setting = "Uniform, n = 150"
A = rbind(A1, A2, A3, A4)

#A$p[A$p == 500] = "p = 500"
#A$p[A$p == 1000] = "p = 1000"
A$p = as.factor(A$p)
A$method = factor(A$method, levels = c("MTT", "QC", "CJ"))
A$p = factor(A$p, levels = c("500", "1000", "1500"))

p3 = A %>% 
  ggplot(aes(y = Power, x = Strength, col = method)) + 
  geom_line(aes(linetype = p)) + geom_point(aes(shape = method)) + 
  facet_grid(.~setting, 
             labeller = labeller(.cols = label_value)) + 
  geom_hline(yintercept = 0.05, colour = "black", linetype = "dashed") + 
  #ggtitle("Uniform distribution") + 
  xlab("Correlation (signal strength)") + 
  scale_x_continuous(breaks = c(0.09, 0.12, 0.15, 0.18, 0.21)) + 
  theme_few() + theme(legend.position = "bottom")  
  #theme(plot.title = element_text(hjust = 0.5, size = 12))

# Figure 3
p3