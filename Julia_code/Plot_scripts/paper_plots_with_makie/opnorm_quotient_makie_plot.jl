using LinearAlgebra, CairoMakie, PolynomialBases, CSV, DelimitedFiles, LaTeXStrings
include("../../functions.jl")
include("../../analysis_tools.jl")

background_norms = readdlm(joinpath(@__DIR__,"./GLL_opnorm_quotient_data/background_norms_data.txt"),'\t', Float64,'\n')
dod_standard_norms = readdlm(joinpath(@__DIR__,"./GLL_opnorm_quotient_data/dod_standard_norms_data.txt"),'\t', Float64,'\n')
dod_optimized_norms = readdlm(joinpath(@__DIR__,"./GLL_opnorm_quotient_data/dod_optimized_norms_data.txt"),'\t', Float64,'\n')


bad_fig_log = Figure(fontsize = 28)
ax_bad = Axis(bad_fig_log[1,1], xlabel = L"p", ylabel = L"Quotient of operator norms$$", xlabelsize = 34, ylabelsize = 34, yscale = log10)
optim_fig = Figure(fontsize = 28)
ax_optim = Axis(optim_fig[1,1], xlabel = L"p", ylabel = L"Quotient of operator norms$$", xlabelsize = 34, ylabelsize = 34)
lines!(ax_bad, range(0, length(dod_standard_norms)-1), vec(dod_standard_norms./background_norms), linewidth = 3.5)
lines!(ax_optim, range(0, length(dod_optimized_norms)-1), vec(dod_optimized_norms./background_norms), linewidth = 3.5)
save(joinpath(@__DIR__,"./opnorm_quotient_bad_log.pdf"), bad_fig_log)
save(joinpath(@__DIR__,"./opnorm_quotient_optim.pdf"), optim_fig)