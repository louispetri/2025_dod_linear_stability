using LinearAlgebra, CairoMakie, PolynomialBases, CSV, DelimitedFiles, LaTeXStrings
include("../../functions.jl")
include("../../analysis_tools.jl")

################################################################################################################################
################################ CODE NOCH NICHT IN PRAXIS GETESTET; KANN FEHLER ENTHALTEN#####################################'
################################################################################################################################
for (basis, lambda_cs, tableind) in 
    [(LobattoLegendre, [1 0.8767 0.5398 0.3203], 1), (GaussLegendre, [1 0.7891 0.4448 0.2781], 2)]
    if basis == LobattoLegendre
        basissstr = "Lobatto"
        shortcut = L"GLL$$"
    else
        basissstr = "Gauss"
        shortcut = L"GL$$"
    end
    for order in [1, 2, 3, 4]
        if basis == LobattoLegendre
            basissstr = "Lobatto"
            titlestring = L"GLL, $p = %$(order-1)$ (order %$(order))"
        else
            basissstr = "Gauss"
            titlestring = L"GL, $p = %$(order-1)$ (order %$(order))"
        end
        f = Figure(fontsize = 24)
        ax = Axis(f[1,1], xlabel = L"Cut-cell factor $\alpha$", ylabel = L"Sharp CFL condition$$", xlabelsize = 34, ylabelsize = 34,
            title = titlestring, titlesize = 36)
        if order in (3,4)
            ylims!(low = 0.0)
        end
        input_timedep = readdlm(joinpath(@__DIR__,"./CFL_lts_data/CFL_lts_timedep_data_$(basissstr)_order=$(order).txt"),'\t', Float64,'\n')
        input_optim = readdlm(joinpath(@__DIR__,"./CFL_lts_data/CFL_lts_optim_data_$(basissstr)_order=$(order).txt"),'\t', Float64,'\n')
        input_background = readdlm(joinpath(@__DIR__,"./CFL_lts_data/CFL_lts_background_data_$(basissstr)_order=$(order).txt"),'\t', Float64,'\n')
        # determine, which dimension is the one distinguising the alphas and the CFL values (because of inconsistent saving)
        if size(input_timedep)[1] < size(input_timedep)[2]
            alpha_res_timedep = input_timedep[1,:]
            CFL_res_timedep = input_timedep[2,:]
        elseif size(input_timedep)[1] > size(input_timedep)[2]
            alpha_res_timedep = input_timedep[:, 1]
            CFL_res_timedep = input_timedep[:, 2]
        else 
            println("CHECK INPUT DATA TIMEDEP")
        end
        if size(input_optim)[1] < size(input_optim)[2]
            alpha_res_optim = input_optim[1,:]
            CFL_res_optim = input_optim[2,:]
        elseif size(input_optim)[1] > size(input_optim)[2]
            alpha_res_optim = input_optim[:, 1]
            CFL_res_optim = input_optim[:,2]
        else 
            println("CHECK INPUT DATA OPTIM")
        end
        if size(input_background)[1] < size(input_background)[2]
            alpha_res_background = input_background[1,:]
            CFL_res_background = input_background[2,:]
        elseif size(input_background)[1] > size(input_background)[2]
            alpha_res_background = input_background[:,1]
            CFL_res_background = input_background[:,2]
        else 
            println("CHECK INPUT DATA BACKGROUND")
        end
        lines!(ax, vec(alpha_res_background), vec(CFL_res_background), label = L"Background without cut cell$", linewidth = 3.5)
        lines!(ax, vec(alpha_res_timedep), vec(CFL_res_timedep), label = L"DoD, $\lambda_c = \lambda_c(\Delta t)$", linestyle = :dash, linewidth = 3.5)
        lines!(ax, vec(alpha_res_optim), vec(CFL_res_optim), label = L"DoD, optimized $\lambda_c$", linestyle = :dashdotdot, linewidth = 3.5)
        #axislegend(position = :lb)
        save(joinpath(@__DIR__,"./CFL_lts_$(basissstr)_order=$(order).pdf"), f)
        g = Figure(fontsize = 20)
        gleg = Legend(g,ax, framevisible = true, labelsize = 17)
        gleg.orientation = :horizontal
        gleg.nbanks = 1
        g[1,1] = gleg
        save(joinpath(@__DIR__,"./legend_CFL_lts.pdf"), g)
    end
end