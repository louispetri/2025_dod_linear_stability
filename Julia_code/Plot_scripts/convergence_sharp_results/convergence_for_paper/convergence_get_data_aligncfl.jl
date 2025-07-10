# uses the table with sharp CFl conditions created in convergence_script.jl to probe the error/convergence error 
# (CFL values are hardcoded inside the function analysis_tools.jl/convergence_error_comparison)
using LinearAlgebra, PolynomialBases, LaTeXStrings, DelimitedFiles, SummationByPartsOperators
include("../../../functions.jl")
include("../../../analysis_tools.jl")


# calculate the CFL values and write in textfile
for basis in [LobattoLegendre, GaussLegendre]
    for order in [1,2,3,4]
        @time steps, errors = convergence_error_comparison(basis, order; exprange = [6,12], sec_factor = 0.95, alignCFL = true, Tmax = 1.0)
        println("done with:$(basis)_order=$(order)")
        # just first row of steps needed, because they are the same for every type of simulation
        data = hcat(steps[:, 1], errors)
			if basis == LobattoLegendre
				basissstr = "Lobatto"
			else
				basissstr = "Gauss"
			end	
        writedlm(joinpath(@__DIR__,"./Sec_fac_0.95_aligncfl/$(basissstr)_order=$(order)_data.txt"), data)
    end
end