# Make a macro that will create matrix R and vector r for restricting the VAR
# The user will write beta[1,2,1] - beta[2,1,3] = 10 (beta[from, to, lag]) and through 
# regular expressions this will get parsed into a vector of the matrix
# macro restrict(ex) gamaa restricts the coefficients

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

# vars, lags, and ntyp have to be defined
# vars = e.vars
# ntyp = e.ntyp
# lags = e.lags

# @restrictions begin
#     gamaa(3)=1
# end

# @restrictions begin
#     gamaa(4) + betaa(1,2,3) + betaa(3,2,1) = 10
# 	betaa(3,4,2) - betaa(3,1,3) = 5
# end


# Heh, hard
# function makeR(R::Matrix{Float64})
# 	base = eye(size(R)[1])

# 	if numberInColumn==1
# 		base

# 	end
# end

# function testRestrictions(R::Matrix{Float64})
# 	@assume rank(R)==apply(max, size(mat)) "The constraints are not independent. One can be expressed as linear combination of other."
# 	# Alright, this is making me nuts :-D
# 	@assume all(mapslices((x) -> sum(map((y) -> y==0.0, x)), getR(mat), [2]) .>= apply(min, size(R))) "Your constants yield the model overidentified. For example, beta(1,2,3) + beta(2,3,1) = 0 and beta(1,2,3) + beta(2,2,1) = 0 imply only equality of beta(2,2,1) - beta(2,3,1) = 0 and does not depend on beta(1,2,3)."
# end