using LinearAlgebra, CairoMakie, PolynomialBases, CSV, DelimitedFiles, LaTeXStrings
include("../../functions.jl")
include("../../analysis_tools.jl")

max_degree = 8;
store_lambda_cs = zeros(max_degree + 1);

c_factors = range(0.001, 0.499, length = 100)
N_spatial = 50
background_norms = zeros(max_degree + 1)
dod_standard_norms = zeros(max_degree + 1)
dod_optimized_norms = zeros(max_degree + 1)
for (ideg, deg) in enumerate(range(0,max_degree))
    for nodes in [LobattoLegendre]
        # Background Operatornrom
        println("Starting degree $(deg)")
        problem = setup_problem_eq("sin", 0, 1, N_spatial, Tmax = 5.0, a = 1.0, CFL = 0.8, bcs = "periodic")
        RHS_mat, problem = DGsemidiscretization_DoD(problem, deg, nodes, "Upwind",  do_stabilize = false)
        background_norms[ideg] = calc_op_norm(problem, RHS_mat)
        # Operatornorm for lambda_c = 1
        dod_standard_norms[ideg] = maximum(op_norm_vs_cut_fac(deg, nodes, c_factors; cellnumber = N_spatial, with_cut_cell = true, do_stabilize = true, fix_eta = true, c = 1))
        # Operatornrm for optimized lambda_c
        (best_norm, lambda_c) = minimize_operatornorm(nodes, deg; cellnumber = N_spatial, a = 1, outflow_cut = true, verbose = false)
        dod_optimized_norms[ideg] = best_norm
        store_lambda_cs .= lambda_c
    end
end
writedlm(joinpath(@__DIR__,"./GLL_opnorm_quotient_data/background_norms_data.txt"), background_norms)
writedlm(joinpath(@__DIR__,"./GLL_opnorm_quotient_data/dod_standard_norms_data.txt"), dod_standard_norms)
writedlm(joinpath(@__DIR__,"./GLL_opnorm_quotient_data/dod_optimized_norms_data.txt"), dod_optimized_norms)
writedlm(joinpath(@__DIR__,"./GLL_opnorm_quotient_data/dod_lambda_data.txt"), store_lambda_cs)