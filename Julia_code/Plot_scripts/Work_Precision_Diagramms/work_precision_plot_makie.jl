using CairoMakie, LaTeXStrings, DelimitedFiles


for basis in [LobattoLegendre, GaussLegendre]
    for order in [1,2,3,4]
        if basis == LobattoLegendre
            basissstr = "Lobatto"
            titlestring = L"GLL, $p = %$(order-1)$ (order %$(order))"
        else
            basissstr = "Gauss"
            titlestring = L"GL, $p = %$(order-1)$ (order %$(order))"
        end
        steps_and_error = readdlm(joinpath(@__DIR__,"./work_precision_data/work_precision_data_$(basissstr)_order=$(order).txt"),'\t', Float64,'\n')
        steps = steps_and_error[1, :]
        stepnumber = length(steps) ÷ 3
        error = steps_and_error[2, :]
        f = Figure(fontsize = 28)
        ax = Axis(f[1,1], xlabel = "time steps", ylabel = "error", yscale = log10, xscale = log10, xlabelsize = 34, ylabelsize = 34,
        title = titlestring,  titlesize = 36)
        lines!(ax, steps[1:stepnumber], error[1:stepnumber], label = L"Background without cut cell$$", linewidth = 3.5)
        lines!(ax, steps[stepnumber+1:2*stepnumber], error[stepnumber+1:2*stepnumber], label = L"DoD, $\lambda_c = \lambda_c(\Delta t)$", linestyle = :dashdotdot, linewidth = 3.5)
        lines!(ax, steps[2*stepnumber+1:3*stepnumber], error[2*stepnumber+1:3*stepnumber], label = L"DoD, optimized $\lambda_c$", linestyle = :dash, linewidth = 3.5)
        #axislegend(ax)
        save(joinpath(@__DIR__,"./work_precision_$(basissstr)_order=$(order).pdf"), f)

        g = Figure(fontsize = 20)
        gleg = Legend(g,ax, framevisible = true, labelsize = 17)
        gleg.orientation = :horizontal
        gleg.nbanks = 1
        g[1,1] = gleg
        save(joinpath(@__DIR__,"./legend_work_precision.pdf"), g)
    end
end
