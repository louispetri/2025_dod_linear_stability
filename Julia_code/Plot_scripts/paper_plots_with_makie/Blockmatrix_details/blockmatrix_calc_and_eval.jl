using LinearAlgebra, CairoMakie, PolynomialBases, LaTeXStrings
include("../../../analysis_tools.jl")

function calc_part_opnorm(setup, RHS_mat, splitRHS; split_ind = 0, ind_col = 0, ind_row = 0)
    Sc = splitRHS["Sc"]
    if split_ind != 10
        col_range = ind_col[split_ind]
        row_range = ind_row[split_ind]
        RHS_part = RHS_mat[row_range, col_range]
    else
        col_range = ind_col[2]
        row_range = ind_row[2]
        RHS_part = splitRHS["DM_J1_just_vol_cm1"][row_range, col_range]
    end
    #Sc[row_range, row_range]
    RHS_part = RHS_mat[row_range, col_range]
    M_global = setup["M"][row_range, row_range]
    dim_M = length(M_global[1,:])
    sqrtM = sqrt.(M_global)
    invsqrtM = zeros(dim_M, dim_M)
    for ism in 1:dim_M
        for jsm in 1:dim_M
            if M_global[ism,jsm] != 0
                invsqrtM[ism,jsm] = 1 ./(sqrt.(M_global[ism,jsm]))
            end
        end
    end
    return opnorm(sqrtM * RHS_part * invsqrtM)
end

# GL nodes, deg = 4
#for (ic, c) in enumerate([1, 0.19529])
# GL nodes, deg = 5
for (deg,optim_c) in [(4,0.19529), (5,0.14927)]
    for (ic, c) in enumerate([1, optim_c])
        cut_cell_number = 3
        cut_cell_sizes = range(0.01, 0.99, length = 200)
        nodes = GaussLegendre
        fix_eta = true
        CFL = 0.2
        do_stabilize = true

        ind_col = 0
        ind_row = 0
        value = zeros(10, length(cut_cell_sizes))
        full_opnorm = zeros(length(cut_cell_sizes))
        for (icut, cut_cell_size) in enumerate(cut_cell_sizes)
            k = deg+1
            # indexes of the influence range: The six kxk block matrices, that influence the cells (c-1, c, c+1).
            ind_infl_range_col = (cut_cell_number-2)*k+1:(cut_cell_number)*k
            ind_infl_range_row = (cut_cell_number-2)*k+1:(cut_cell_number+1)*k

            # indexes of the block matrixes that appear within the estimate
            ind_L_cm1L_col = (cut_cell_number-3)*k+1:(cut_cell_number-2)*k
            ind_L_cm1_col = (cut_cell_number-2)*k+1:(cut_cell_number-1)*k
            ind_L_cm1R_col = (cut_cell_number-1)*k+1:(cut_cell_number)*k
            ind_L_cL_col = (cut_cell_number-2)*k+1:(cut_cell_number-1)*k
            ind_L_c_col = (cut_cell_number-1)*k+1:(cut_cell_number)*k
            ind_L_cp1LL_col = (cut_cell_number-2)*k+1:(cut_cell_number-1)*k
            ind_L_cp1L_col = (cut_cell_number-1)*k+1:(cut_cell_number)*k
            ind_L_cp1_col = (cut_cell_number)*k+1:(cut_cell_number+1)*k
            ind_L_cp2L_col = (cut_cell_number)*k+1:(cut_cell_number+1)*k
            ind_L_cm1L_row = (cut_cell_number-2)*k+1:(cut_cell_number-1)*k
            ind_L_cm1_row = (cut_cell_number-2)*k+1:(cut_cell_number-1)*k
            ind_L_cm1R_row = (cut_cell_number-2)*k+1:(cut_cell_number-1)*k
            ind_L_cL_row = (cut_cell_number-1)*k+1:(cut_cell_number)*k
            ind_L_c_row = (cut_cell_number-1)*k+1:(cut_cell_number)*k
            ind_L_cp1LL_row = (cut_cell_number)*k+1:(cut_cell_number+1)*k
            ind_L_cp1L_row = (cut_cell_number)*k+1:(cut_cell_number+1)*k
            ind_L_cp1_row = (cut_cell_number)*k+1:(cut_cell_number+1)*k
            ind_L_cp2L_row = (cut_cell_number+1)*k+1:(cut_cell_number+2)*k
            
            ind_col = [ind_L_cm1L_col, ind_L_cm1_col, ind_L_cm1R_col, ind_L_cL_col, ind_L_c_col, ind_L_cp1LL_col, ind_L_cp1L_col, ind_L_cp1_col, ind_L_cp2L_col]
            ind_row = [ind_L_cm1L_row, ind_L_cm1_row, ind_L_cm1R_row, ind_L_cL_row, ind_L_c_row, ind_L_cp1LL_row, ind_L_cp1L_row, ind_L_cp1_row, ind_L_cp2L_row]
            
            problem = setup_problem_eq("sin", 0, 1, 50, Tmax = 5, a = 1, CFL = CFL);
            problem = include_cut_cell(problem, cut_cell_size, cut_cell_number);
            RHS_mat, setup, splitRHS = DGsemidiscretization_DoD(problem, deg, nodes, "Upwind", do_stabilize = do_stabilize, c=c, fix_eta = fix_eta);
            for i in 1:10
                value[i, icut] = calc_part_opnorm(setup, RHS_mat, splitRHS, split_ind = i, ind_col = ind_col, ind_row = ind_row)
            end        
        end
        full_opnorm = op_norm_vs_cut_fac(deg, nodes, cut_cell_sizes; with_cut_cell = true, do_stabilize = do_stabilize, cutnumber = cut_cell_number, fix_eta = fix_eta, outflow_cut = true, CFL = CFL, c = c);

        labellist = [L"L_{(c-1)L}", L"L_{(c-1)}", L"L_{(c-1)R}", L"L_{cL}", L"L_{c}", L"L_{(c+1)LL}", L"L_{(c+1)L}"
                        , L"L_{(c+1)}", L"L_{(c+2)L}", L"(\star)"]

        f = Figure(fontsize = 23)
        if ic == 1
            ax = Axis(f[1,1], xlabel = L"Cut-cell factor $\alpha$", ylabel = L"Operator norm$$", title = L"p=%$(deg), \;\lambda_c = %$(c)", xlabelsize = 34, ylabelsize = 28, titlesize = 36)
        else
            ax = Axis(f[1,1], xlabel = L"Cut-cell factor $\alpha$", title =  L"p=%$(deg), \;\lambda_c = %$(c)", xlabelsize = 34, titlesize = 36)
        end
        for i in 1:10
            if i != 10 && i != 9
                lines!(ax, cut_cell_sizes, abs.(value[i,:]), label = labellist[i], linewidth = 3)
            elseif i == 9
                lines!(ax, cut_cell_sizes, abs.(value[i,:]), label = labellist[i], linewidth = 3, color = :purple)
            else # i == 10
                lines!(ax, cut_cell_sizes, abs.(value[i,:]), label = labellist[i], linestyle = :dash, linewidth = 3, color = :green2)
            end
        end
        lines!(ax, cut_cell_sizes, full_opnorm, label = L"L", color = :firebrick2, linestyle = :dash, linewidth = 3)
        xlims!(low = 0, high = 0.49)
        if ic == 1
            ylims!(low = 000, high =  40000)
        end
        save(joinpath(@__DIR__,"./partial_opnorm_c=$(c)_GL_deg=$(deg).pdf"), f)

        g = Figure(fontsize = 23)
        gleg = Legend(g[1,1],ax, framevisible = true, position = (1.5, 1.5))
        gleg.orientation = :horizontal
        gleg.nbanks = 2
        save(joinpath(@__DIR__,"./partial_opnorm_legend.pdf"), g)
    end
end

#gleg = Legend(f[1,3],ax, framevisible = true, fontsize = 23)

# (*) = L"S_{c-1}M^{-1}\eta_c\mathbb{I}^TD^TM\mathbb{I}"
