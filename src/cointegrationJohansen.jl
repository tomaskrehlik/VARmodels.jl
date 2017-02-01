function criticalValues(ecdet, typ)
	if ecdet=="none"
		cvals = [6.5, 12.91, 18.9, 24.78, 30.84, 36.25, 42.06, 48.43, 54.01, 59, 65.07, 8.18, 14.9, 21.07, 27.14, 33.32, 39.43, 44.91, 51.07, 57, 62.42, 68.27, 11.65, 19.19, 25.75, 32.14, 38.78, 44.59, 51.3, 57.07, 63.37, 68.61, 74.36, 6.5, 15.66, 28.71, 45.23, 66.49, 85.18, 118.99, 151.38, 186.54, 226.34, 269.53, 8.18, 17.95, 31.52, 48.28, 70.6, 90.39, 124.25, 157.11, 192.84, 232.49, 277.39, 11.65, 23.52, 37.22, 55.43, 78.87, 104.2, 136.06, 168.92, 204.79, 246.27, 292.65]		
	elseif ecdet=="const"
		cvals = [7.52, 13.75, 19.77, 25.56, 31.66, 37.45, 43.25, 48.91, 54.35, 60.25, 66.02, 9.24, 15.67, 22, 28.14, 34.4, 40.3, 46.45, 52, 57.42, 63.57, 69.74, 12.97, 20.2, 26.81, 33.24, 39.79, 46.82, 51.91, 57.95, 63.71, 69.94, 76.63, 7.52, 17.85, 32, 49.65, 71.86, 97.18, 126.58, 159.48, 196.37, 236.54, 282.45, 9.24, 19.96, 34.91, 53.12, 76.07, 102.14, 131.7, 165.58, 202.92, 244.15, 291.4, 12.97, 24.6, 41.07, 60.16, 84.45, 111.01, 143.09, 177.2, 215.74, 257.68, 307.64]
	elseif ecdet=="trend"
		cvals = [10.49, 16.85, 23.11, 29.12, 34.75, 40.91, 46.32, 52.16, 57.87, 63.18, 69.26, 12.25, 18.96, 25.54, 31.46, 37.52, 43.97, 49.42, 55.5, 61.29, 66.23, 72.72, 16.26, 23.65, 30.34, 36.65, 42.36, 49.51, 54.71, 62.46, 67.88, 73.73, 79.23, 10.49, 22.76, 39.06, 59.14, 83.2, 110.42, 141.01, 176.67, 215.17, 256.72, 303.13, 12.25, 25.32, 42.44, 62.99, 87.31, 114.9, 146.76, 182.82, 222.21, 263.42, 310.81, 16.26, 30.45, 48.45, 70.05, 96.58, 124.75, 158.49, 196.08, 234.41, 279.07, 327.45]
	end
	cvals = reshape(cvals, 11, 3, 2)
	if typ=="trace"
		return convert(Matrix, cvals[:,:,1])
	elseif typ=="eigen"
		return convert(Matrix, cvals[:,:,2])
	end
end

function selector(i, n)
	v = zeros(n)
	v[i] = 1
	return v
end

function cointegrationJohansen(data, ecdet, K, spec; season = "no", dumvar = "no")
	@assert K >= 2 "K must be at least 2."
	N, P = size(data)
	season == "no" ? s = 0 : s = season - 1
	@assert N * P > P + s * P + K * P^2 + P * (P + 1)/2 "Insufficient degrees of freedom."

	
	@assert ecdet in ["none"; "const"; "trend"] "Unsupported structure of equation, ecdet must be one of: none, const, trend"
	@assert spec in ["longrun"; "transitory"] "Specification of equation not supported, spec must be one of longrun, transitory"

	nanindexes = !mapslices(x->any(isnan(x)), data, [2])

	if season != "no"
		seasondums = (hcat([selector(mod(x,season)+1,season) for x=1:N]...)')[:,2:season]
		seasondums = seasondums[nanindexes, :]
	end

	if dumvar != "no"
		@assert size(dumvar)[1]==N "Unequal row length between dummy variables and data matrix."
	end

	if ecdet == "trend"
		trend = collect(1:N)
		trend = trend[nanindexes, :]
	end

	data = data[nanindexes, :]
    N, P = size(data)
    Z = lagData(vcat(mapslices(x->diff(x), data, [1]), zeros(P)'), K, N)
    Z0 = Z[:, 1:P]

    
    spec == "longrun" ? Lnotation = K : Lnotation = 1

    if (ecdet == "none") && (spec == "longrun")
    	ZK = data[1:(N - K), :]
    	Z1 = Z[:, (P+1):(size(Z)[2])]
        Z1 = hcat([1.0 for i in 1:size(Z1)[1]], Z1)
        idx = collect(0:(P - 1))
        model = "with linear trend"
    end
    if (ecdet == "none") && (spec == "transitory")
		ZK = data[1:(N-1), :][K:(N - 1), :]
		Z1 = Z[:, (P+1):(size(Z)[2])]
        Z1 = hcat([1.0 for i in 1:size(Z1)[1]], Z1)
        idx = collect(0:(P - 1))
        model = "with linear trend"
    end
    if (ecdet == "const") && (spec == "longrun")
    	ZK = hcat(data[1:(N - K), :], [1 for i in 1:(N-K)])
    	Z1 = Z[:, (P+1):(size(Z)[2])]
    	P = P + 1
        idx = collect(0:(P - 2))
        model = "without linear trend and constant in cointegration"
    end
    if (ecdet == "const") && (spec == "transitory")
    	ZK = hcat(data[1:(N-1), :], [1 for i in 1:(N-1)])[K:(N - 1), :]
    	Z1 = Z[:, (P+1):(size(Z)[2])]
    	P = P + 1
        idx = collect(0:(P - 2))
        model = "without linear trend and constant in cointegration"
    end
    if (ecdet == "trend") && (spec == "longrun")
		ZK = hcat(data[1:(N - K), :], trend[1:(N - K)])
		Z1 = Z[:, (P+1):(size(Z)[2])]
        Z1 = hcat([1.0 for i in 1:size(Z1)[1]], Z1)
        P = P + 1
        idx = collect(0:(P - 2))
        model = "with linear trend in cointegration"
    end
    if (ecdet == "trend") && (spec == "transitory")
    	ZK = hcat(data[1:(N-1), :], trend[1:(N-1)])[K:(N - 1), :]
    	Z1 = Z[:, (P+1):(size(Z)[2])]
        Z1 = hcat([1.0 for i in 1:size(Z1)[1]], Z1)
        P = P + 1
        idx = collect(0:(P - 2))
        model = "with linear trend in cointegration"
    end
    N = size(Z0)[1]

    if season != "no"
    	if ecdet == "const"
    		Z1 = hcat(seasondums[(K+1):size(seasondums)[1],:], Z1)
    	else
    		Z1 = hcat(Z1[:, 1], seasondums[(K+1):size(seasondums)[1],:], Z1[:, 2:size(Z1)[2]])
    	end
    end

    if dumvar != "no"
    	if ecdet == "const"
    		Z1 = hcat(dumvar[(K+1):size(dumvar)[1],:], Z1)
    	else
    		Z1 = hcat(Z1[:, 1], dumvar[(K+1):size(dumvar)[1],:], Z1[:, 2:size(Z1)[2]])
    	end
    end



    M00 = (Z0'*Z0)/N
    M11 = (Z1'*Z1)/N
    MKK = (ZK'*ZK)/N
    M01 = (Z0'*Z1)/N
    M0K = (Z0'*ZK)/N
    MK0 = (ZK'*Z0)/N
    M10 = (Z1'*Z0)/N
    M1K = (Z1'*ZK)/N
    MK1 = (ZK'*Z1)/N
    M11inv = inv(M11)
    R0 = Z0 - (M01 * M11inv * Z1')'
    RK = ZK - (MK1 * M11inv * Z1')'
    S00 = M00 - M01 * M11inv * M10
    S0K = M0K - M01 * M11inv * M1K
    SK0 = MK0 - MK1 * M11inv * M10
    SKK = MKK - MK1 * M11inv * M1K
    Ctemp = cholfact(SKK, :U, Val{true})
    pivot = sortperm(Ctemp.piv)
    C = (convert(Matrix,Ctemp[:U])[:, pivot])'
    Cinv = inv(C)
    S00inv = inv(S00)
    λ, e = eig(Cinv * SK0 * S00inv * S0K * Cinv')
    V = Cinv' * e
    Vorg = V
    V = hcat([V[:, x]./V[1, x] for x in 1:P]...)
    W = S0K * V * inv(V' * SKK * V)
    Π = S0K * inv(SKK)
    Δ = S00 - S0K * V * inv(V' * SKK * V) * V' * SK0
    Γ = M01 * M11inv - Π * MK1 * M11inv

    return (data, Z0, Z1, ZK, dumvar, λ, Vorg, V, W, Π, Δ, Γ, R0, RK, spec, season)
end

function johansenTeststat(cointEst, typ)
	@assert typ in ["eigen"; "trace"] "Unsupported type, typ must be one of: eigen, trace"
	P = size(cointEst.ZK)[2]
	N = size(cointEst.data)[1]

	if (ecdet == "none") && (spec == "longrun")
        idx = collect(0:(P - 1))
    end
    if (ecdet == "none") && (spec == "transitory")
        idx = collect(0:(P - 1))
    end
    if (ecdet == "const") && (spec == "longrun")
        idx = collect(0:(P - 2))
    end
    if (ecdet == "const") && (spec == "transitory")
        idx = collect(0:(P - 2))
    end
    if (ecdet == "trend") && (spec == "longrun")
        idx = collect(0:(P - 2))
    end
    if (ecdet == "trend") && (spec == "transitory")
        idx = collect(0:(P - 2))
    end

	if P > 11
		warning("Too many variables, critical values cannot be computed.")
		cvals = []
	else
		cvals = criticalValues(cointEst.ecdet, typ)
	end	

	if typ == "trace"
    	teststat = reverse(map(x -> -N * sum(log(1-sort(cointEst.λ, rev = true)[x:P])), idx + 1))
    elseif typ == "eigen"
    	teststat = reverse( (-N*log(1-sort(cointEst.λ, rev = true)))[idx + 1])
    end
    return (teststat, cvals)
end

type CointegrationEstimate
	data::Matrix{Float64}
	Z0::Matrix{Float64}
	Z1::Matrix{Float64}
	ZK::Matrix{Float64}
    dumvar
    λ::Vector{Float64}
    Vorg::Matrix{Float64}
    V::Matrix{Float64}
    W::Matrix{Float64}
    PI::Matrix{Float64}
    DELTA::Matrix{Float64}
    GAMMA::Matrix{Float64}
    R0::Matrix{Float64}
    RK::Matrix{Float64}
    spec::AbstractString
    season
    ecdet::AbstractString

	function CointegrationEstimate(data::Matrix{Float64}, ecdet::AbstractString, K::Int, spec::AbstractString; season = "no", dumvar = "no")
		data, Z0, Z1, ZK, dumvar, λ, Vorg, V, W, PI, DELTA, GAMMA, R0, RK, spec, season = cointegrationJohansen(data, ecdet, K, spec; season = season, dumvar = dumvar)
		new(data, Z0, Z1, ZK, dumvar, λ, Vorg, V, W, PI, DELTA, GAMMA, R0, RK, spec, season, ecdet)
	end
end

function getVARrepresentation(cointObject, r::Int)
	@assert typeof(cointObject)==CointegrationEstimate "The first argument should be cointegration estimate."
	@assert r < size(cointObject.data)[2] "The cointegration rank---the second argument---can be no bigger than number of variables in the system."

	N, vars = size(cointObject.data)


	etc = cointObject.ZK * cointObject.V[:, 1:r]
	coeffs = hcat(etc, cointObject.Z1) \ cointObject.Z0

	Π = cointObject.W[:, 1:r] * cointObject.V[:,1:r]'

	if cointObject.ecdet == "const"
		detcoeffs = Π[:, vars + 1]
		Π = Π[:, 1:vars]
		rhs = hcat([1 for i in 1:size(cointObject.ZK)[1]], cointObject.Z1)
	elseif cointObject.ecdet == "none"
		detcoeffs = coeffs[r+1,:]'
		rhs = cointObject.Z1
	elseif cointObject.ecdet == "trend"
		detcoeffs = hcat(coeffs[r+1,:]',Π[:,vars + 1])
		Π = Π[:, 1:vars]
		rhs = hcat([1 for i in 1:size(cointObject.ZK)[1]], cointObject.ZK[:, vars + 1], cointObject.Z1[:, 2:size(cointObject.Z1)[2]])
	end

	cointObject.ecdet == "const" ? a = 0 : a = 1
	if cointObject.season != "no"
		s = cointObject.season - 1
		coefs_from_coeffs = coeffs[(r+a) + 1:s, :]
		detcoeffs = vcat(detcoeffs, coefs_from_coeffs)
	else
		s = 0
	end
	if cointObject.dumvar != "no"
		d = size(cointObject.dumvar)[2]
		coefs_from_coeffs = coeffs[(r+a+s) + 1:d, :] 
		detcoeffs = vcat(detcoeffs, coefs_from_coeffs)
	else
		d = 0
	end

	detcoeffs = detcoeffs'
	
	lags = N-size(cointObject.ZK)[1]
	Γ = coeffs[((-(lags-1)*vars+1):0)+size(coeffs)[1],:]'
	A = zeros(vars, vars, lags)

	i = hcat([vars*i + 1 for i=0:(lags-2)], [vars + vars*i for i=0:(lags-2)])
	if cointObject.spec == "transitory"
		A[:,:,1] = Γ[:, i[1,1]:i[1,2]] + Π + eye(vars)
		for j=2:(lags-1)
			A[:,:,j] = Γ[:, i[j,1]:i[j,2]] - Γ[:,i[j - 1,1]:i[j - 1,2]]
		end
		A[:,:,lags] = -Γ[:, i[lags-1,1]:i[lags-1, 2]]
	elseif cointObject.spec == "longrun"
		A[:,:,1] = Γ[:, i[1,1]:i[1,2]] + eye(vars)
		for j=2:(lags-1)
			A[:,:,j] = Γ[:, i[j,1]:i[j,2]] - Γ[:,i[j - 1,1]:i[j - 1,2]]
		end
		A[:,:,lags] = Π-Γ[:, i[lags-1,1]:i[lags-1, 2]]
	end

	determin_data = rhs[:,1:size(detcoeffs)[1]]
	data = lagData(vcat(cointObject.data, zeros(vars)'), lags, N)
	Y = (cointObject.data)[(lags+1):size(cointObject.data)[1],:]
	C = vcat(detcoeffs, vcat([A[:,:,i]' for i=1:size(A)[3]]...))
	X = hcat(determin_data, data)
	
	resids = Y - X * C

	Σ = cov(resids)*(N-1)/(N-lags*vars-size(detcoeffs)[1]) 
	return (lags, size(detcoeffs)[1], vars, N, X, Y, C, Σ)
end

type varRepresentationVECM <: VARRepresentation
	lags::Int
	typ::AbstractString
	ntyp::Int
	vars::Int
	obs::Int
	X::Matrix{Float64}
	Y::Matrix{Float64}
	C::Matrix{Float64}
	Σ::Matrix{Float64}
	r::Int
	HPsi::Int
	HPhi::Int
	Phi::Array{Float64, 3}
	Psi::Array{Float64, 3}
	seriesNames::Vector{AbstractString}

	function varRepresentationVECM(cointegrationEstimate, r::Int)
		(lags, ntyp, vars, obs, X, Y, C, Σ) = getVARrepresentation(cointegrationEstimate, r::Int)
		new(lags, "varRepresentation", ntyp, vars, obs, X, Y, C, Σ, r, 0, 0, zeros(vars,vars,1), zeros(vars,vars,1), ["" for i=1:vars])
	end
end
