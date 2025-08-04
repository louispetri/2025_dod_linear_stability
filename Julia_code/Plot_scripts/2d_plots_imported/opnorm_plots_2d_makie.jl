using LinearAlgebra, CairoMakie, PolynomialBases, DelimitedFiles, LaTeXStrings
include("../../functions.jl")
include("../../analysis_tools.jl")



# Operator norm plots for 2d
background_data = readdlm(joinpath(@__DIR__,"./Plotdaten/operator-norm-background"),' ', Float64,'\n')
for (deg, yranges, lambdas, leg_pos) in [(2, [200, 1050], (0.1925, 0.5, 1.0), :lt), (3, [200, 4500], (0.0875, 0.3, 1.0), :rt), (4, [000, 14000], (0.046, 0.2, 1.0), :rt)]
    f = Figure(fontsize = 24)
    ax = Axis(f[1,1], xlabel = L"Cut-cell factor $\alpha$", ylabel = L"Operator norm$$", title = L"p=%$deg", xlabelsize = 34, ylabelsize = 34, titlesize = 36, ytickformat = values -> ["$(Int(value))" for value in values])
    xlims!(ax, low = 0, high = 0.42)
    ylims!(ax, low = yranges[1], high = yranges[2])
    linestyles = [:solid, :dash, :dashdotdot]
    background_alphas = zeros(1)
    for i in 1:3
        import_data = readdlm(joinpath(@__DIR__,"./Plotdaten/operator-norm-order-$(deg)-tau-$(lambdas[i])"),' ', Float64,'\n')
        line = lines!(ax, import_data[:, 1], import_data[:, 2], label = L"$\lambda_c = %$(lambdas[i])$", linestyle = linestyles[i], linewidth = 3.5)
        if i == 1
            background_alphas = import_data[:, 1]
        end
    end
    line = lines!(ax, background_alphas, ones(size(background_alphas))*background_data[deg-1, 2], label = L"$$Background", linestyle = :dashdotdot, linewidth = 3.5)
    axislegend(ax, position = leg_pos)
    save(joinpath(@__DIR__,"./opnorm_2d_p=$(deg)_makie.pdf"), f)
end