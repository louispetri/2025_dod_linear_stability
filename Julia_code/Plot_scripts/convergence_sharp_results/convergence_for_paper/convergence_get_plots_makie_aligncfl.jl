# uses the table with sharp CFl conditions created in convergence_script.jl to probe the error/convergence error 
# (CFL values are hardcoded inside the function analysis_tools.jl/convergence_error_comparison)
using LinearAlgebra, CairoMakie, PolynomialBases, LaTeXStrings, DelimitedFiles, SummationByPartsOperators
include("../../../functions.jl")
include("../../../analysis_tools.jl")

# plot the data

for basis in [LobattoLegendre, GaussLegendre]
    for order in [1,2,3,4]
        if basis == LobattoLegendre
            basissstr = "Lobatto"
            titlestring = L"GLL, $p = %$(order-1)$ (order %$(order))"
        else
            basissstr = "Gauss"
            titlestring = L"GL, $p = %$(order-1)$ (order %$(order))"
        end
        B = readdlm(joinpath(@__DIR__,"./Sec_fac_0.95_aligncfl/$(basissstr)_order=$(order)_data.txt"),'\t', Float64,'\n')

        f = Figure(fontsize = 28)
        ax = Axis(f[1,1], xlabel = L"1/\Delta x", ylabel = L"error$$", xlabelsize = 34, ylabelsize = 34,
            title = titlestring, titlesize = 36, xscale = log10, yscale = log10)
        lines!(ax, B[:,1], B[:,2], label = L"Background without cut cell$$", linewidth = 3.5)
        lines!(ax, B[:,1], B[:,3], label = L"DoD, $\lambda_c = \lambda_c(\Delta t)$", linestyle = :dash, linewidth = 3.5)
        lines!(ax, B[:,1], B[:,4], label = L"DoD, optimized $\lambda_c$", linestyle = :dashdotdot, linewidth = 3.5)
        save(joinpath(@__DIR__,"./Sec_fac_0.95_aligncfl/$(basissstr)_order=$(order)_plot.pdf"), f)

        g = Figure(fontsize = 20)
        gleg = Legend(g,ax, framevisible = true, labelsize = 17)
        gleg.orientation = :horizontal
        gleg.nbanks = 1
        g[1,1] = gleg
        save(joinpath(@__DIR__,"./Sec_fac_0.95_aligncfl/legend_convergence.pdf"), g)
    end
end