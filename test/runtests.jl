using VARmodels
using Base.Test

# write your own tests here
i = rand(1:2000)
run(`Rscript varest.R $i`)
C = readcsv("coefs.csv")
S = readcsv("sigma.csv")
resids = readcsv("residuals.csv")

data = readcsv("data.csv")
data = log(data)
est = varEstimate(data[(1+i):(500+i),:], 2, "Const")
residuals = (est.Y-est.X*est.C)
@test all(abs(C[2:size(C)[1],:]-est.C).<1e-10)
@test all(abs(residuals-resids[2:size(resids)[1],:]).<1e-10)
@test all(abs(S[2:size(S)[1],:]-est.Σ).<1e-10)
run(`rm coefs.csv`)
run(`rm residuals.csv`)
run(`rm sigma.csv`)