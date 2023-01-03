# Optimal-Test-Cov

**functions.R** includes the main R functions for the proposed MT test (*MTT*), the maximal type test (*CJ*) proposed in Cai and Jiang (2011) and the sum-of-square type test (*QC*) proposed in Qiu and Chen (2012). The functions *MTT* and *CJ* are written in R, the function *QC* needs to call C functions. 

**bandcov.R** includes the R functions to call C for the L2 test statistic in Qiu and Chen (2012).

**bandcov.c** includes the C functions for the L2 test statistic in Qiu and Chen (2012).

**Compile C functions:**
