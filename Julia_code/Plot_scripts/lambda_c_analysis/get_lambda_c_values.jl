# Calculate the optimized values for lambda_c up to an arbitrary degree for both types of nodes and store them in a text file
using LinearAlgebra, PolynomialBases, DelimitedFiles, LaTeXStrings
include("../../functions.jl")
include("../../analysis_tools.jl")

# Choose highest number of polynomial degree
deg_max = 12
degs = 0:deg_max
cellnumber = 50
# Store in row (one, three) the (lambda_c, opnorm) - values for LobattoLegendre and in row (two, four) the (lambda_c, opnorm) - values for GaussLegendre nodes
c_values = zeros(length(degs),2)
opnorm_values = zeros(length(degs),2)
for basis in [LobattoLegendre, GaussLegendre]
    for (ideg, deg) in enumerate(degs)
        opnormval, cval = minimize_operatornorm(basis, deg, cellnumber = cellnumber, a = 1, outflow_cut = true)
        if basis == LobattoLegendre
            c_values[ideg, 1] = cval
            opnorm_values[ideg, 1] = opnormval
        else 
            c_values[ideg, 2] = cval
            opnorm_values[ideg, 2] = opnormval
        end
        println("------------------STAGE $(deg) COMPLETED--------------------")
    end
end
data = hcat(c_values, opnorm_values)
writedlm(joinpath(@__DIR__,"./opt_lambda_c&opnorm_xn=$(cellnumber).txt"), data)