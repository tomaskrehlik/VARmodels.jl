# All the VAR.jl, SVAR.jl, SVEC.jl (will) have the following components
# - estimateXXXX function
# - XXXXEstimate type
# - λ function meaning the log-likelihood
# - Phi to get the Wold representation
# - Psi to get the orthogonalized Wold representation
# - show function that gives the output

# Additional machinery such as tests should be properly referenced to other files

# In the long run, the following should be supported:
# - Analytical plots using the PyPlot library
# - Automatic selection using information criteria (or at least giving table with the numbers)
# - Macro language for the restricted VAR using the GLS.


module VARmodels
	using Formatting
	using Distributions

	include("dataManipulations.jl")
	include("VAR.jl")
	include("restrictions.jl")
	include("bindings.jl")
	include("fevd.jl")
	include("plots.jl")

	# Main VAR machinery
	export 
		varEstimate,
		Psi,
		restrictVAR,
		restrictVAR2,
		getR,
		Phi

	# FEVD functions
	export 
		fevd,
		genFEVD,
		fftFEVD,
		fftGenFEVD

	# Simulation functions
	export
		simulateVAR

	# Plotting
	export
		plotIRF


end # module
