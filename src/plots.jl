function plotIRF(estimate, H)
	irfs = Psi(estimate, H)
	n, n, H = size(irfs)
	fig = figure("pyplot_subplot_mixed",figsize=(10,10))
	for j = 1:n
		for i = 1:n
			subplot(n,n,(i-1)*n+j)
			grid("on")
			plot(vec(irfs[i,j,:]))
		end
	end
	fig[:canvas][:draw]()
	fig.close()
end