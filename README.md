# Optimal-Test-Cov

**functions.R** includes the main R functions for the proposed MT test (*MTT*), the maximal type test (*CJ*) proposed in Cai and Jiang (2011) and the sum-of-square type test (*QC*) proposed in Qiu and Chen (2012). The functions *MTT* and *CJ* are written in R, the function *QC* needs to call C functions. 

**bandcov.R** includes the R functions to call C for the L2 test statistic in Qiu and Chen (2012).

**bandcov.c** includes the C functions to compute the L2 test statistic in Qiu and Chen (2012).

**Compile C functions:** Run the following line in R to generate ".so" file (for Mac and server) or to generate ".dll" file (for Windows). 

system("R CMD SHLIB -o bandcov.so bandcov.c")

For Windows systems, please replace the line *dyn.load("bandcov.so")* in the function *QC* in the file **functions.R** by *dyn.load("bandcov.dll")*.

**size-normal.R**, **size-normal-p1500.R**, **size-uniform.R** and **size-uniform-p1500.R** are the simulation codes for the sizes of the three tests under normal and uniform distributed data. Those four codes are used to make Table 1 in the paper.

**power-normal-n100.R**, **power-normal-n100-p1500.R**, **power-normal-n150.R**, **power-normal-n150-p1500.R**, **power-uniform-n100.R**, **power-uniform-n100-p1500.R**, **power-uniform-n150.R** and **power-uniform-n150-p1500.R** are the simulation codes for the power of the three tests under normal and uniform distributed data, and sparse signal settings.

**plot-new.R** is the R code to make Figures 2 and 3 in the paper for the powers of the three tests.
