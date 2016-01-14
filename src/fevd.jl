function fevd(estimate::varEstimate, H)
	Psi(estimate, H)
	FEVD = [sum(estimate.Psi[i,j,:].^2) for i=1:estimate.vars, j=1:estimate.vars]
	FEVD = hcat([FEVD[:,i]/sum(FEVD[:,i]) for i=1:estimate.vars]...)
	return FEVD
end

function genFEVD(estimate::varEstimate, H::Int, corzer::Bool = false)
	Phi(estimate, H)

	A = estimate.Phi
	Σ = deepcopy(estimate.Σ)
	K = estimate.vars

	den = diag(mapreduce((i)->(A[:,:,i]'*Σ*A[:,:,i]), +, collect(1:H)))

	if corzer
		Σ = diagm(diag(Σ))
		num = mapreduce((i)->(A[:,:,i]'*Σ).^2, +, collect(1:H))
	else
		num = mapreduce((i)->(A[:,:,i]'*Σ).^2, +, collect(1:H))	
	end

	θ = [num[i,j]/(den[i]*sqrt(Σ[j,j])) for i=1:K, j=1:K]
	for i=1:K
		θ[i,:] = θ[i,:]/sum(θ[i,:])
	end
	θ = convert(Array{Float64}, θ)

	return θ
end

# Fourier decomposition of the FEVDs

function fftFEVD(estimate::varEstimate, H)
	k = estimate.vars

	Phi(estimate, H)

	P = chol(estimate.Σ)

	ft = mapslices(fft, estimate.Phi, [3])
	decomp = [(abs(ft[:,i,z]'*P'[:,j]).^2)[1] / H for i=1:k, j = 1:k, z = 1:H]
	denom = mapslices(sum, decomp, [3, 2])

	return  mapslices((x) -> x./denom[:,:,1], decomp, [1,2])
end

function fftGenFEVD(estimate::varEstimate, H, nocorr = false)
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
		for j=1:k
			for h=1:H
				decompNew[i,j,h] = (abs(ft[:,i,h]'*Σ'[:,j]).^2)[1] / H
			end
		end
	end
	# decompNew = [(abs(ft[:,i,h]'*Σ'[:,j]).^2)[1] / H for i=1:k, j = 1:k, h = 1:H]
	decomp = zeros(k,k,H)
	for i=1:k
		for j=1:k
			for h=1:H
				decomp[i,j,h] = (abs(ft[:,i,h]'*P'[:,j]).^2)[1] / H
			end
		end
	end
	# decomp = [(abs(ft[:,i,h]'*P'[:,j]).^2)[1] / H for i=1:k, j = 1:k, h = 1:H]
	denom = vec(mapslices(sum, decomp, [3, 2]))
	
	θ = [decompNew[i,j,h]/(denom[i]*sqrt(Σ[j,j])) for i=1:k, j = 1:k, h = 1:H]
	# display(θ)
	# div = vec(mapslices(sum, θ, [2,3]))

	for i=1:k
		s = sum(θ[i,:,:])
		for j=1:k
			for h=1:H
				θ[i,j,h] = θ[i,j,h]/s
			end
		end
	end

	# θ = [θ[i,j,h] / sum(θ[i,:,:]) for i=1:k, j = 1:k, h = 1:H]

	return θ
end