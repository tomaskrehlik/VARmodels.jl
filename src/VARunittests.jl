using Distributions

# Unit-tests for the VAR estimation

lags = 3
series = 6
length = 1000

#
# # Function for simulating Vector autoregressive process with given coefficient matrix
#

# This function simulates Vector autoregressive process taking those values as input:
# - length: length of the resulting time series
# - typ: type of the process, whether it should or should not have a constant, etc.
# - coefficients: matrix of coefficients with dimensions (series, series*lags) with ordering in each line as x_lag1, y_lag1, x_lag2, y_lag2, etc...
# - other_c: matrix of coefficients corresponding to potential constants and trends with size (series, ...) where ... is 0, 1, or 2, according to type

function VARsim(length, typ, coefficients, other_c)
	# Should be enough
	burnout = 200

	a, b = size(coefficients)
	lags = convert(Int, max(a,b)/min(a,b))
	series = min(a,b)

	typ=="None" && (other=fill(0.0, length + burnout))
	typ=="Const" && (other=fill(1.0, length + burnout))
	typ=="Trend" && (other=[1:(length+burnout)])
	typ=="Const and trend" && (other=[fill(1.0, length + burnout) [1:(length+burnout)]])

	y = fill(0.0, length + burnout, series)
	eps = rand(Normal(0, 1/100), (length + burnout, series))

	for i=(lags+1):(size(y)[1])
			y[i,:] = mapreduce((x) -> y[i-x,:]*coefficients[(1:series) + (x-1)*series, :], +, 1:lags) + eps[i,:] + other[i,:]*other_c'
	end

	y = y[(burnout+1):end,:]
	
	return y
end

# 
# # Function simulating random VAR of given specification
# 
# This function is suited for performing a unit-test of the VAR system.
# Arguments:
# - lags: number of lags in the system
# - series: number of series in the system
# - length: length of the resulting time series
# - typ: whether the system should not have constant, trend or have any combination of those.

function VARdatarandsim(lags, series, length, typ)
	typ=="None" && (other_c=fill(0.0, series))
	typ=="Const" && (other_c=rand(Uniform(), series))
	typ=="Trend" && (other_c=rand(Uniform(), series))
	typ=="Const and trend" && (other_c=[rand(Uniform(), series) rand(Uniform(), series)])

	# Set them so that it does not blow up!
	coefficients = rand(Uniform(0, 1/(lags*series)), (series, lags*series))'
	y = VARsim(length, typ, coefficients, other_c)
	return (y, [other_c', coefficients])
end



apply(vcat,[coeficients[(1:series) + (x-1)*series, :] for x=[lags:-1:1]]) - varEstimate(y, lags, "None").C

mapreduce((x)->VARdifcoef(lags, series, length, burnout), +, 1:60)/60