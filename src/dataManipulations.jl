#
###### Data transformations
#

# Lag data takes in a matrix. If we had columns a b c, then we will have a matrix with
# a(-lags) b(-lags) c(-lags) a(-lags + 1) b(-lags + 1) c(-lags + 1) ... ... ... a(-1) b(-1) c(-1).
# Maybe the columns should be reordered to have the closest lags at the front? changed!

function lagData(data, lags::Int64, obs::Int64)
	VARdata = data[lags:(obs-1),:]
	for i=(lags-1):-1:1
		VARdata = [VARdata data[i:(obs-lags+i-1),:]]
	end
	return VARdata
end

# dataMatrix is a function that creates a data matrix that goes into the estimation
# args:
#	- data: matrix of numbers
# 	- lags: number of lags in the estimation
# 	- typ: ASCIIString describing the type of VAR to be fitted, can be either of the following
# 			- None, do not include neither constant nor trend
# 			- Const, do include a constant
# 			- Trend, do include a trend
# 			- Const and trend, do include both constant and trend
# Returns a data matrix as from lagData with first columns corresponding to type.

function dataMatrix(data, lags, typ)
	obs = size(data)[1]
	data = lagData(data, lags, obs)

	typ=="None" && return(data)
	typ=="Const" && return([ones(obs-lags) data])
	typ=="Const and trend" && return([ones(obs-lags) 1:(obs-lags) data])
	typ=="Trend" && return([1:(obs-lags) data])
end

#
######
#
