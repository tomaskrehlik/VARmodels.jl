library(vars)

i <- as.integer(commandArgs(trailingOnly = TRUE)[1])

data <- read.csv("data.csv", header = F)
data <- log(data)

est <- VAR(data[1:500+i, ], 2, "const")
C <- sapply(est$varresult, function(i) coefficients(i))
C <- C[c(nrow(C),1:(nrow(C)-1)),]
write.csv(file = "coefs.csv", C, row.names = F)
write.csv(file = "sigma.csv", summary(est)$cov, row.names = F)
write.csv(file = "residuals.csv", residuals(est), row.names = F)