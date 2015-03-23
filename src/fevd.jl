function fevd(estimate::varEstimate, H)
	Psi(estimate, H)
	FEVD = [sum(estimate.Psi[i,j,:].^2) for i=1:estimate.vars, j=1:estimate.vars]
	FEVD = apply(hcat,[FEVD[:,i]/sum(FEVD[:,i]) for i=1:estimate.vars])
	return FEVD
end

function genFEVD(estimate::varEstimate, H::Int, corzer::Bool = false)
	Phi(estimate, H)

	A = estimate.Phi
	Σ = deepcopy(estimate.Σ)
	K = estimate.vars

	den = diag(mapreduce((i)->(A[:,:,i]'*Σ*A[:,:,i]), +, [1:H]))

	if corzer
		Σ = diagm(diag(Σ))
		num = mapreduce((i)->(A[:,:,i]'*Σ).^2, +, [1:H])
	else
		num = mapreduce((i)->(A[:,:,i]'*Σ).^2, +, [1:H])	
	end

	θ = [num[i,j]/(den[i]*Σ[i,i]) for i=1:K, j=1:K]
	for i=1:K
		θ[:,i] = θ[:,i]/sum(θ[:,i])
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

	P = chol(Σ)

	ft = mapslices(fft, estimate.Phi, [3])

	decompNew = [(abs(ft[:,i,z]'*Σ'[:,j]).^2)[1] / H for i=1:k, j = 1:k, z = 1:H]
	decomp = [(abs(ft[:,i,z]'*P'[:,j]).^2)[1] / H for i=1:k, j = 1:k, z = 1:H]
	denom = mapslices(sum, decomp, [3, 2])

	θ = zeros(k,k,H)

	for i = 1:k
		for j = 1:k
			for h = 1:H
				θ[i,j,h] = decompNew[i,j,h]/(denom[i,1,1]*Σ[i,i])
			end
		end
	end

	div = mapslices(sum, θ, [1,3])

	for i = 1:k
		for j = 1:k
			for h = 1:H
				θ[i,j,h] = θ[i,j,h] / div[1,j,1]
			end
		end
	end

	return θ
end