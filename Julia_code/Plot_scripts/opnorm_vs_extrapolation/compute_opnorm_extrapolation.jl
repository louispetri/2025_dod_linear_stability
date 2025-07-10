using PolynomialBases, LinearAlgebra, CairoMakie, LaTeXStrings

f = Figure(fontsize = 18)
g = Figure(fontsize = 18)
ax1GLL = Axis(f[1,1], xlabel = L"\alpha", ylabel = L"\Vert\mathbb{I}\Vert^{1/p}", title = L"GL$",
        xlabelsize = 28, ylabelsize = 28, titlesize = 28)
ax1GL = Axis(f[2,1], xlabel = L"\alpha", ylabel = L"\Vert\mathbb{I}\Vert^{1/p}", title = L"GLL$",
xlabelsize = 28, ylabelsize = 28, titlesize = 28)
alphas = range(0,0.49, length = 100)
opnorms = zeros(length(alphas))
for p in range(2,15)
    for nodes in [LobattoLegendre, GaussLegendre]
        basis  = nodes(p)
        sqrt_mass = diagm(sqrt.(basis.weights))
        invsqrt_mass = diagm(sqrt.(1 ./basis.weights))
        global opnorms = [opnorm(sqrt_mass*interpolation_matrix(1 .+ alpha .+ alpha * basis.nodes, basis)*invsqrt_mass) for alpha in alphas]
        if nodes == LobattoLegendre
            lines!(ax1GLL, alphas, opnorms.^(1/(p)), label = "$(p)")
        else
            lines!(ax1GL, alphas, opnorms.^(1/(p)), label = "$(p)")
        end
    end
end
f[1:2, 2] = Legend(f, ax1GLL, L"orders$", framevisible = true)
#axislegend(ax1GL, position = :lt)
#axislegend(ax1GLL, position = :lt)

axsqrt = Axis(g[2,1], xlabel = L"p", ylabel = L"\Vert\mathbb{I}\Vert^{1/p}",
xlabelsize = 28, ylabelsize = 28, titlesize = 28)
axlog = Axis(g[1,1], xlabel = L"p", ylabel = L"\Vert\mathbb{I}\Vert", yscale = log10,
xlabelsize = 28, ylabelsize = 28, titlesize = 28)
xlims!(axsqrt, low = 1, high = 30)
xlims!(axlog, low = 1, high = 30)
ylims!(axsqrt, low = 1, high = 4)
ps = range(0,30)
alpha = 0.4
opnorms2 = zeros(length(ps))
for nodes in [LobattoLegendre, GaussLegendre]
    for (pi, p) in enumerate(ps)
        basis  = nodes(p)
        sqrt_mass = diagm(sqrt.(basis.weights))
        invsqrt_mass = diagm(sqrt.(1 ./basis.weights))
        opnorms2[pi] = opnorm(sqrt_mass*interpolation_matrix(1 .+ alpha .+ alpha * basis.nodes, basis)*invsqrt_mass)
    end
    op_exp = 1 ./ps
    if nodes == LobattoLegendre
        lines!(axsqrt, ps, opnorms2 .^op_exp, label = L"GLL$", linewidth = 1.5; color = :blue)
        lines!(axlog, ps, opnorms2 , label = L"GLL$", linewidth = 1.5; color = :blue)
    else
        lines!(axsqrt, ps, opnorms2 .^op_exp, label = L"GL$", linestyle = :dash, linewidth = 2.5; color = :tomato)
        lines!(axlog, ps, opnorms2 , label = L"GL$", linestyle = :dash, linewidth = 2.5; color = :tomato)
    end
end
#glegend = Legend(g, axsqrt, "node types", framevisible = true)
#glegend.nbanks = 1
#g[1:2, 2] = glegend
axislegend(axsqrt, position = :rb)
axislegend(axlog, position = :rb)



save(joinpath(@__DIR__,"./extrapol_alpha_vs_opnorm.pdf"), f)
save(joinpath(@__DIR__,"./extrapol_p_vs_opnorm.pdf"), g)

# auch daran denken, noch plots wie in arbitrary order details einzubringen. D.h. p-opnorm plots