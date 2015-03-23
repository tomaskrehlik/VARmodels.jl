module VARmodels
	using Formatting

	include("dataManipulations.jl")
	include("VAR.jl")
	include("restrictions.jl")
	include("bindings.jl")
	include("fevd.jl")

	export 
		varEstimate,
		fevd,
		genFEVD,
		fftFEVD,
		fftGenFEVD,
		Psi,
		restrictVAR,
		restrictVAR2,
		getR,
		Phi

end # module
