# 
# # Example with VAR as in Sims (1980)
# 

using VARmodels, TimeSeries

data = readdlm("/Users/tomaskrehlik/.julia/v0.3/VARmodels/src/simsData.txt", ';')

d = [Date(1959,1,1):Month(3):Date(2014,7,1)]
names = convert(Vector{String},vec(data[1,2:end]))

data = TimeArray(d, convert(Matrix{Float64},data[46:end,2:end]), names)

s = varEstimate(data["wages","imports","cpi"], 3, "Const")




Series: ....
Lags: ....
Type: ....

No restrictions imposed. Estimated by OLS.

