using LaTeXStrings, DelimitedFiles, PolynomialBases
function work_precision_diagramm(basis, order; sec_factor = 1.0, alignCFL = false, relative_path_folder = "../convergence_sharp_results/No_securty_factor")
    # using the same CFL values as for the convergence tests and also Tmax = 5, a=1. (All values hardcoded in convergence_error_comparison and its called functions)
    CFL_backgr = [1 1 0.44999 0.74672 ; 1 0.33481 0.2098 0.449994]
    CFL_DoD_classic = [1 .789668 0.164004 0.0976313 ; 1 0.323104 0.119105 0.0815261]
    CFL_opt = [1 0.78966 0.27625 0.358243 ; 1 0.32700 0.152292 0.26454]
    if alignCFL == true
        CFL_min = zeros(length(CFL_opt[1,:]))
        for i in range(1,length(CFL_opt[1,:]))
            CFL_min[i] = minimum([CFL_backgr[i], CFL_DoD_classic[i], CFL_opt[i]])
            CFL_backgr[i] = CFL_min[i]
            CFL_DoD_classic[i] = CFL_min[i]
            CFL_opt[i] = CFL_min[i]
        end
    end
    # import data: B[:,1] = spatialsteps, B[:,2] = background method, B[:,3] = classic DoD, B[:,4] = optimized DoD
    if basis == LobattoLegendre
        basissstr = "Lobatto"
    else
        basissstr = "Gauss"
    end	
    if basis == LobattoLegendre
        nodes_index = 1
    elseif basis == GaussLegendre
        nodes_index = 2
    end
    B = readdlm(joinpath(@__DIR__, relative_path_folder * "/$(basissstr)_order=$(order)_data.txt"),'\t', Float64,'\n')
    # determine done timesteps:
    C_backgr = CFL_backgr[nodes_index, order]
    C_DoD_classic = CFL_DoD_classic[nodes_index, order]
    C_opt = CFL_opt[nodes_index, order]
    tsteps_backgr = [round(5*xsteps/(sec_factor*C_backgr))-1 for xsteps in B[:,1]]
    tsteps_DoD_classic = [round(5*xsteps/(sec_factor*C_DoD_classic))-1 for xsteps in B[:,1]]
    tsteps_opt = [round(5*xsteps/(sec_factor*C_opt))-1 for xsteps in B[:,1]]

    tsteps = hcat(tsteps_backgr, tsteps_DoD_classic, tsteps_opt)
    return tsteps, B[:,2:4]
end