function fevd(estimate::VARRepresentation, H)
	Psi(estimate, H)
	FEVD = [sum(estimate.Psi[i,j,:].^2) for i=1:estimate.vars, j=1:estimate.vars]
	FEVD = hcat([FEVD[:,i]/sum(FEVD[:,i]) for i=1:estimate.vars]...)
	return FEVD
end

function genFEVD(estimate::VARRepresentation, H::Int; nocorr = false)
	Phi(estimate, H)

	A = estimate.Phi
	Σ = deepcopy(estimate.Σ)
	K = estimate.vars

	den = zeros(size(A)[1])
	num = zeros(A[:,:,1])
	
	if nocorr
		Σ = diagm(diag(deepcopy(estimate.Σ)))
	end

	for i=1:H
		t = A[:,:,i]
		z = t'*Σ
		den += diag(z*t)
		num += z.^2
	end
	
	θ = [num[i,j]/(den[i]*Σ[j,j]) for i=1:K, j=1:K]
	for i=1:K
		θ[i,:] = θ[i,:]/sum(θ[i,:])
	end
	θ = convert(Array{Float64}, θ)

	return θ
end

# Fourier decomposition of the FEVDs

function fftFEVD(estimate::VARRepresentation, H; nocorr = false, range = nothing)
	k = estimate.vars

	Phi(estimate, H)

	if nocorr
		Σ = diagm(diag(deepcopy(estimate.Σ)))
	else
		Σ = estimate.Σ
	end

	P = chol(Σ)

	if !(typeof(range) <: Tuple)
		range[1] = 1
		range[2] = H
	end

	ft = mapslices(fft, estimate.Phi, [3])
	decomp = [(abs(ft[:,i,z]'*P'[:,j]).^2)[1] / H for i=1:k, j = 1:k, z = 1:H]
	denom = mapslices(sum, decomp[:,:,range[1]:range[2]], [3, 2])

	return  mapslices((x) -> x./denom[:,:,1], decomp, [1,2])
end

function fftGenFEVD(estimate::VARRepresentation, H; nocorr = false, range = nothing)
	k = estimate.vars

	Phi(estimate, H)

	if nocorr
		Σ = diagm(diag(estimate.Σ))
	else
		Σ = estimate.Σ
	end

	P = convert(Array{Float64,2}, chol(Σ))
	ft = mapslices(fft, estimate.Phi, [3])
	decompNew = zeros(k,k,H)
	for i=1:k
		temp = ft[:,i,:]
		for j=1:k
			temp2 = Σ[j,:]
			for h=1:H
				decompNew[i,j,h] = (abs(temp[:,h]'*temp2).^2)[1] / H
			end
		end
	end

	decomp = zeros(k,k,H)

	for i=1:k
		for j=1:k
			for h=1:H
				decomp[i,j,h] = (abs(ft[:,i,h]'*P'[:,j]).^2)[1] / H
			end
		end
	end

	if !(typeof(range) <: Tuple)
		range = [1,H]
	end		

	denom = vec(mapslices(sum, decomp[:,:,range[1]:range[2]], [3, 2]))

	θ = [decompNew[i,j,h]/(denom[i]*Σ[j,j]) for i=1:k, j = 1:k, h = 1:H]

	for i=1:k
		s = sum(θ[i,:,range[1]:range[2]])
		for j=1:k
			for h=1:H
				θ[i,j,h] = θ[i,j,h]/s
			end
		end
	end

	return θ
end
