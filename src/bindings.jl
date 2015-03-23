import TimeSeries: TimeArray

function varEstimate(data::TimeArray{Float64}, lags::Int, typ::String)
	t = varEstimate(data.values, lags, typ)
	t.seriesNames = data.colnames
	return t
end