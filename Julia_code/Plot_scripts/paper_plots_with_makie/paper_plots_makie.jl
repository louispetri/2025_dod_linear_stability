using LinearAlgebra, CairoMakie, PolynomialBases, CSV, DelimitedFiles, LaTeXStrings
include("../../functions.jl")
include("../../analysis_tools.jl")



# Operator norm plots for old DoD with bad operatornrom
for (deg, yranges) in [(2, [000, 5000]), (3, [000, 5000]), (4, [000, 8000])]
    f = Figure(fontsize = 24)
    ax = Axis(f[1,1], xlabel = L"Cut-cell factor $\alpha$", ylabel = L"Operator norm$$", title = L"p=%$deg", xlabelsize = 34, ylabelsize = 34, titlesize = 36)
    xlims!(ax, low = 0, high = 0.49)
    ylims!(ax, low = yranges[1], high = yranges[2])
    CFL = 0.5
    nodes = LobattoLegendre
    c_factors = range(0.01, 0.99, length = 200)
    cutnum = 3
    outflow_cut = true
    
    norms = op_norm_vs_cut_fac(deg, nodes, c_factors; with_cut_cell = false, do_stabilize = false, CFL = CFL, outflow_cut = outflow_cut)
    background = lines!(ax, c_factors, norms, label = L"Background without cut cell$$", linewidth = 3.5)
    
    norms = op_norm_vs_cut_fac(deg, nodes, c_factors; with_cut_cell = true, do_stabilize = false, cutnumber = cutnum, CFL = CFL, outflow_cut = outflow_cut)
    without_stabil = lines!(ax, c_factors, norms, label = L"Background with cut cell$$", linewidth = 3.5)
    
    norms = op_norm_vs_cut_fac(deg, nodes, c_factors; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, CFL = CFL, outflow_cut = outflow_cut)
    stabil_c1 = lines!(ax, c_factors, norms, label = L"DoD, $\lambda_c = 1$", linestyle = :dash, linewidth = 3.5)
    norms = op_norm_vs_cut_fac(deg, nodes, c_factors; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, CFL = CFL, outflow_cut = outflow_cut, c= 0.5)
    stabil_c05 = lines!(ax, c_factors, norms, label = L"DoD, $\lambda_c = 0.5$", linestyle = :dash, color = :mediumorchid1, linewidth = 3.5)
    norms = op_norm_vs_cut_fac(deg, nodes, c_factors; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, CFL = CFL, outflow_cut = outflow_cut, c= 2.0)
    stabil_c2 = lines!(ax, c_factors, norms, label = L"DoD, $\lambda_c = 2$", linestyle = :dash, color = :firebrick2, linewidth = 3.5)
    save(joinpath(@__DIR__,"./p=$(deg)_Lobatto_bad_operatornorm_plot.pdf"), f)
    g = Figure(fontsize = 20)
    gleg = Legend(g,ax, framevisible = true)
    gleg.orientation = :horizontal
    gleg.nbanks = 2
    g[1,1] = gleg
    save(joinpath(@__DIR__,"./legend_Lobatto_bad_operatornorm_plot.pdf"), g)
end



# Operatr norm plots for DoD with optimized operator norm
for (deg, yranges, lambda) in [(2, [000, 5000], 0.5399), (3, [000, 5000], 0.3213), (4, [000, 8000], 0.2230)]
    f = Figure(fontsize = 24)
    ax = Axis(f[1,1], xlabel = L"Cut-cell factor $\alpha$", ylabel = L"Operator norm$$",title = L"$p=%$(deg),\;$ optimized $ \lambda_c = %$(lambda)$" , xlabelsize = 34, ylabelsize = 34, titlesize = 36)
    xlims!(ax, low = 0, high = 0.49)
    ylims!(ax, low = yranges[1], high = yranges[2])
    CFL = 0.5
    nodes = LobattoLegendre
    c_factors = range(0.01, 0.99, length = 200)
    cutnum = 3
    outflow_cut = true
    
    norms = op_norm_vs_cut_fac(deg, nodes, c_factors; with_cut_cell = false, do_stabilize = false, CFL = CFL, outflow_cut = outflow_cut)
    background = lines!(ax, c_factors, norms, label = L"Background without cut cell$$", linewidth = 3.5)
    
    norms = op_norm_vs_cut_fac(deg, nodes, c_factors; with_cut_cell = true, do_stabilize = false, cutnumber = cutnum, CFL = CFL, outflow_cut = outflow_cut)
    without_stabil = lines!(ax, c_factors, norms, label = L"Background with cut cell$$", linewidth = 3.5)

    norms = op_norm_vs_cut_fac(deg, nodes, c_factors; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, CFL = CFL, outflow_cut = outflow_cut)
    stabil_c1 = lines!(ax, c_factors, norms, label = L"DoD, $\lambda_c = 1$", linestyle = :dash, linewidth = 3.5)
    
    norms = op_norm_vs_cut_fac(deg, nodes, c_factors; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, CFL = CFL, outflow_cut = outflow_cut, c = lambda)
    stabil_coptim = lines!(ax, c_factors, norms, label = L"DoD, $\lambda_c$ optimized", linestyle = :dash, color = :mediumorchid1, linewidth = 3.5)

    save(joinpath(@__DIR__,"./p=$(deg)_Lobatto_optim_operatornorm_plot.pdf"), f)
    g = Figure(fontsize = 20)
    gleg = Legend(g,ax, framevisible = true)
    gleg.orientation = :horizontal
    gleg.nbanks = 2
    g[1,1] = gleg
    save(joinpath(@__DIR__,"./legend_Lobatto_optim_operatornorm_plot.pdf"), g)
end