using LinearAlgebra, CairoMakie, PolynomialBases, CSV, DelimitedFiles, LaTeXStrings
include("../../functions.jl")
include("../../analysis_tools.jl")


CFL_conditions = zeros(8,3);
for (basis, lambda_cs, tableind) in 
	[(LobattoLegendre, [1 0.8767 0.5398 0.3203], 1), (GaussLegendre, [1 0.7891 0.4448 0.2781], 2)]
    if basis == LobattoLegendre
        basissstr = "Lobatto"
    else
        basissstr = "Gauss"
    end
		for (order, TMM) in [(1, expl_Euler), (2, Heun), (3,SSPRK3), (4, SSPRK10_4)]
			outflow_cut = true
			#plot(xlabel = "cut-cell-factor", ylabel = "sharp CFL condition")
			
			alpha_res, CFL_res = calc_lts_CFL(basis, order - 1, TMM, cellnumber = 50, a=1, background = false, fix_eta = false, outflow_cut = outflow_cut)
            writedlm(joinpath(@__DIR__,"./CFL_lts_data/CFL_lts_timedep_data_$(basissstr)_order=$(order).txt"), [alpha_res, CFL_res])

			#plot!(alpha_res, CFL_res, label = "DoD, c=c(Delta t)")
			CFL_conditions[(tableind-1)*4+order, 1] = minimum(CFL_res)
			
			alpha_res, CFL_res = calc_lts_CFL(basis, order - 1, TMM, cellnumber = 50, a=1, background = false, fix_eta = true, outflow_cut = outflow_cut, c = lambda_cs[order])
			writedlm(joinpath(@__DIR__,"./CFL_lts_data/CFL_lts_optim_data_$(basissstr)_order=$(order).txt"), [alpha_res, CFL_res])
            #plot!(alpha_res, CFL_res, label = "DoD, lambda_c=$(lambda_cs[order])")
			CFL_conditions[(tableind-1)*4+order, 2] = minimum(CFL_res)
			
			alpha_res, CFL_res = calc_lts_CFL(basis, order -1, TMM, cellnumber = 50, a=1, background = true)
			writedlm(joinpath(@__DIR__,"./CFL_lts_data/CFL_lts_background_data_$(basissstr)_order=$(order).txt"), [alpha_res CFL_res])
            #plot!(alpha_res, CFL_res, label = "background", linestyle = :dashdot)
			CFL_conditions[(tableind-1)*4+order, 3] = minimum(CFL_res)			
			#savefig(plot!(), "./CFL_LTS_true_labeled/$(basissstr)_order=$(order).pdf")
            println("Saved Data4")
		end
	end
    writedlm(joinpath(@__DIR__,"./CFL_lts_data/CFL_lts_c(deltat)_data.txt"), [CFL_conditions[:,3] CFL_conditions[:,2] CFL_conditions[:,1]])
	#df = DataFrame(method = ["GLB1", "GLB2", "GLB3", "GLB4", "GLL1", "GLL2", "GLL3", "GLL4"], CFL_backgr = CFL_conditions[:,3], CFL_optim = CFL_conditions[:,2], CFL_timedep = CFL_conditions[:,1])
	#plt = plot(;legend=nothing,xaxis=false,yaxis=false,xticks=false,yticks=false)
	#plt = annotate!(0.0,1.0,(PrettyTables.pretty_table(String,df),10,:left,:top,:black,:Courier))
	#savefig(plt,"./CFL_LTS_true_labeled/table_CFL_cond.pdf")