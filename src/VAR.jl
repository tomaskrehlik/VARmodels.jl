# This script should eventually become a ultimate copy of the all encompasing 
# textbook on multivariate time-series: Lütkepohl (2007). Some of the methods will be added
# from other papers, such as generalised VAR structure from Pesaran, Shin or computing 
# all reorderings in the standard identification structure.

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

#
###### Definition of the specialised type
# 

# The definition of the type rests on several specific neccessities
# There are things that are neccessary for proper definition of the VAR model and as
# such have to be included in the definition of the type. That is for example number
# of lags, the original data, the type.
# Other things can be foregone, but are often desirable to be stored, because they are
# computationally intensive and computing them only once is typically much better than
# recomputing them every time. However, in many applications, one might not even need them.
# 
# It is neccessary to balance the needs and strike the balance between initial speed of
# estimation and future speed. For example, I probably wouldn't like to have bootstrapped
# confidence intervals for some things, but I would surely keep the possibility to have them.
# 
# Moreover, the current estimation technique will probably be replaced by the method from 
# function restrictVAR2, because it is more general.

function estimateVAR(data::Matrix{Float64}, lags::Int, typ::String)
	X = dataMatrix(data, lags, typ)
	Y = data[(1+lags):end,:]
	obs, vars = size(data)
	C = (X'*X) \ X'*Y

	# Get numerical type
	typ=="None" && (ntyp = 0)
	typ=="Const" && (ntyp = 1)
	typ=="Const and trend" && (ntyp = 2)
	typ=="Trend" && (ntyp = 1)

	# Σ = (X*C-Y)'*(X*C-Y) / (obs-lags-lags*vars-1)
	Σ = (X*C-Y)'*(X*C-Y) / (obs-lags)

	return (X, Y, obs, vars, C, Σ, ntyp)
end

type varEstimate
	lags::Int
	typ::ASCIIString
	ntyp::Int
	vars::Int
	obs::Int
	X::Matrix{Float64}
	Y::Matrix{Float64}
	C::Matrix{Float64}
	Σ::Matrix{Float64}
	HPsi::Int
	HPhi::Int
	Phi::Array{Float64, 3}
	Psi::Array{Float64, 3}

	function varEstimate(data::Matrix{Float64},lags::Int,typ::String)
		(X, Y, obs, vars, C, Σ, ntyp) = estimateVAR(data::Matrix{Float64},lags::Int,typ::String)
		new(lags, typ, ntyp, vars, obs, X, Y, C, Σ, 0, 0, zeros(vars,vars,1), zeros(vars,vars,1))
	end

end

# type varRestrictions
# 	R::Matrix{Float64}
# 	r::Vector{Float64}
# 	desc::Vector{ASCIIString}

# 	function varRestrictions(R::Matrix{Float64}, r::Vector{Float64}, desc::Vector{ASCIIString})
# 		@assume size(R)[1] == size(r)
# 		new(R, r, desc)
# 	end
# end

function Clag(est::varEstimate, lag)
	C = est.C
	ntyp = est.ntyp
	(coefs, vars) = size(C)
	lags = (coefs-ntyp)/vars
	lags < lag ? error("The demanded lag is higher than existing") : ""
	return C[1+ntyp+(((lag-1)*vars):((lag)*vars-1)),:]
end

function getR(a::Matrix{Bool})
	b = prod(size(a))
	c = sum(a)
	mat = fill(0.0, b, c)
	a = vec(a)
	d=1
	for i=1:b
		if a[i]
			mat[i,:] = [fill(0.0, d-1), 1, fill(0.0, c-d)]
			d += 1
		else
			mat[i,:] = fill(0.0, c)
		end
	end

	return mat
end

# Make a macro that will create matrix R and vector r for restricting the VAR
# The user will write beta[1,2,1] - beta[2,1,3] = 10 (beta[from, to, lag]) and through 
# regular expressions this will get parsed into a vector of the matrix
# macro restrict(ex)

function β(from, to, lag)
	col = fill(0.0, (vars^2)*lags + vars*ntyp)
	println(1 + vars*ntyp + (lag-1)*(vars^2) + (from-1)*vars + (to-1))
	col[1 + vars*ntyp + (lag-1)*(vars^2) + (from-1)*vars + (to-1)] = 1.0
	return col
end

function betaa(from, to, lag)
	return β(from, to, lag)
end

function γ(eq)
	col = fill(0.0, (vars^2)*lags + vars*ntyp)
	col[eq] = 1.0
	return col
end

function gamaa(eq)
	return γ(eq)
end

macro restrictions(ex)
	ind = [2:2:length(ex.args)]
	R = apply(Expr, prepend!([:($(ex.args[i].args[1])) for i=ind], [:call, :hcat]));
	r = [:($(ex.args[i].args[2].args[2])) for i=ind];
	return :(tuple($R, $r))
end

@restrictions begin
    gamaa(3)=1
end

# Heh, hard
# function makeR(R::Matrix{Float64})
# 	base = eye(size(R)[1])

# 	if numberInColumn==1
# 		base

# 	end
# end

function testRestrictions(R::Matrix{Float64})
	@assume rank(R)==apply(max, size(mat)) "The constraints are not independent. One can be expressed as linear combination of other."
	# Alright, this is making me nuts :-D
	@assume all(mapslices((x) -> sum(map((y) -> y==0.0, x)), getR(mat), [2]) .>= apply(min, size(R))) "Your constants yield the model overidentified. For example, beta(1,2,3) + beta(2,3,1) = 0 and beta(1,2,3) + beta(2,2,1) = 0 imply only equality of beta(2,2,1) - beta(2,3,1) = 0 and does not depend on beta(1,2,3)."
end

function restrictVAR(est::varEstimate, restrictions::Expr, egls::Bool=false)
	# TODO: Add small r from Lutkepohl 5.2
	# The original data are not changed, due to very different structure in the estimation (vectorisation)
	# The matrix R as the argument has to be of the size of the coefficients
	# Y stay the same

	# Naming conventions, this does not eat up memory, it is just renaming for more readability, 
	# it is parsed out
	Z = est.X'
	Y = est.Y
	Σ = est.Σ
	L = est.lags
	K = est.vars
	k = est.obs-est.lags
	R = restrictions.R
	r = restrictions.r

	z = vec(Y) - kron(eye(K), Z')*r

	if egls
		Σ1 = inv(Σ)
		γ = R'*kron(Σ1, Z*Z')*R \ R'*kron(Σ1, Z)*z
	else
		γ = (R'*kron(eye(K), Z*Z')*R) \ R'*kron(eye(K), Z)*z
	end
	β = R*γ + r
	est.C = reshape(β, L*K+1, K)

	est.Σ = (est.X*est.C-est.Y)'*(est.X*est.C-est.Y) / k
	return est
end


function Phi(estimate::varEstimate, H)
	phi = zeros(estimate.vars, estimate.vars, H)
	phi[:,:,1] = eye(estimate.vars)

	for j=2:H
		k = j > estimate.lags ? estimate.lags + 1 : j
		for i=1:(k-1)
			phi[:,:,j] += phi[:,:,j-i]*Clag(estimate, i)
	 	end
	end
	estimate.Phi = phi
	estimate.HPhi = H
	return phi
end

function Psi(estimate::varEstimate, H)
	Phi(estimate, H)
	psi = deepcopy(estimate.Phi)
	P = chol(estimate.Σ)
	for i=1:H
		psi[:,:,i] = P*psi[:,:,i]
	end
	estimate.Psi = psi
	estimate.HPsi = H
	return psi
end

function fevd(estimate::varEstimate, H)
	Psi(estimate, H)
	FEVD = [sum(estimate.Psi[i,j,:].^2) for i=1:estimate.vars, j=1:estimate.vars]
	FEVD = apply(hcat,[FEVD[:,i]/sum(FEVD[:,i]) for i=1:estimate.vars])
	return FEVD
end

