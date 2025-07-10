#############################################################################
############################ Convergence test ###############################
#############################################################################


function test_convergence_background_method(basis, order ; TMM = SSPRK10_4, exprange = [3,8], CFL = 0.9, a=1, Tmax = 5)
    stepsizes = 2 .^range(exprange[1],exprange[2])
    errors = zeros(length(stepsizes))
    for (istep, stepsize) in enumerate(stepsizes)
        problem = setup_problem_eq("sin", 0, 1, stepsize, Tmax = Tmax, a = a, CFL = CFL)
        RHS_mat, problem = DGsemidiscretization_DoD(problem, order, basis, "Upwind", do_stabilize = false)
        solution, u_exact = TMM(problem, RHS_mat)
        #errors[istep] = norm(solution[:, end] - u_exact[:, end])
         
        x_d = problem["x_d"]
        v = problem["v"]
        deg = problem["deg"]
        nodes = problem["nodes"]
        errors[istep] = 0
        for i in 1:problem["cellnumber"]
            basis1 = nodes(deg)
            for j in 1:deg+1
                # Inserting these nodes into the structure basis1, to adapt the field basis1.interpolation_matrix to the used nodes.
                basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
            end
            errors[istep] += integrate((solution[(deg+1)*(i-1)+1:(deg+1)*i, end] - u_exact[(deg+1)*(i-1)+1:(deg+1)*i, end]).^2, basis1.weights.*(v[i+1]-v[i])/2)
        end
        errors[istep] = sqrt(errors[istep])
        
    end
    return stepsizes, errors
end


# Due to unexpected better convergecne rates: Test with more/other nodes (solved now, so this is not needed anymore, but still could be a nice check)
function test_convergence_background_method1(basis, order ; exprange = [3,8], CFL = 0.9, a=1)
    stepsizes = 2 .^range(exprange[1],exprange[2])
    errors = zeros(length(stepsizes))
    solution = 0
    solution_ev = 0
    x_d = 0
    x_d_ev = 0
    for (istep, stepsize) in enumerate(stepsizes)
        problem = setup_problem_eq("sin", 0, 1, stepsize, Tmax = 5, a = a, CFL = CFL)
        RHS_mat, problem = DGsemidiscretization_DoD(problem, order, basis, "Upwind", do_stabilize = false)
        solution, u_exact = SSPRK10_4(problem, RHS_mat)
        #errors[istep] = norm(solution[:, end] - u_exact[:, end])
         
        x_d = problem["x_d"]
        v = problem["v"]
        deg = problem["deg"]
        nodes = problem["nodes"]
        cellnumber = problem["cellnumber"]
        t_d = problem["t_d"]
        u0_eval = problem["u0_eval"]
        shift(i) = (v[i+1]+v[i])/2
        scale(i) = (v[i+1]-v[i])/2
        
        # Interpolate solution at 2p+1 nodes instead p+1 = q  -> _ev = 2*q - 1 = 2p+2 - 1 = 2p+1.
        deg_ev = 2*deg
        x_d_ev = zeros(cellnumber*(deg_ev+1))
        solution_ev = zeros(cellnumber*(deg_ev+1))
        basis_ev = GaussLegendre(deg_ev)
        basis1 = nodes(deg)
        for i in 1:cellnumber
            for j in 1:deg_ev+1
                x_d_ev[(deg_ev+1)*(i-1) + j] = basis_ev.nodes[j]*scale(i) + shift(i)
            end
        end
        u_exact_ev = [u0_eval(to_periodic(x-(a*(t)))) for x in x_d_ev, t in t_d]

        for i in 1:cellnumber
            for j in 1:deg+1
                basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
            end
            transfo = interpolation_matrix(x_d_ev[(deg_ev+1)*(i-1) + 1:(deg_ev+1)*(i)], basis1)
            solution_ev[(deg_ev+1)*(i-1) + 1:(deg_ev+1)*(i)] = transfo * solution[(deg+1)*(i-1)+1:(deg+1)*i]
        end
        
        
        
        errors[istep] = 0
        
        for i in 1:cellnumber
            basis_ev1 = GaussLegendre(deg_ev)
            for j in 1:deg_ev+1
                # Inserting these nodes into the structure basis1, to adapt the field basis1.interpolation_matrix to the used nodes.
                basis_ev1.nodes[j] = x_d_ev[(deg_ev+1)*(i-1) + j] 
            end
            errors[istep] += integrate((solution_ev[(deg_ev+1)*(i-1)+1:(deg_ev+1)*i, end] - u_exact_ev[(deg_ev+1)*(i-1)+1:(deg_ev+1)*i, end]).^2, basis_ev1.weights.*(v[i+1]-v[i])/2)
        end
        errors[istep] = sqrt(errors[istep])
        
        
    end
    return stepsizes, errors, solution, solution_ev, x_d, x_d_ev
end

function test_convergence_dod_method(basis, order ; exprange = [3,8], TMM = SSPRK10_4, CFL = 0.9, a=1, outflow_cut = true, fix_eta = false, c=1, cut_cell_factor = 0.01, Tmax = 5)
    stepsizes = 2 .^range(exprange[1],exprange[2])
    errors = zeros(length(stepsizes))
    for (istep, stepsize) in enumerate(stepsizes)
        problem = setup_problem_eq("sin", 0, 1, stepsize, Tmax = Tmax, a = a, CFL = CFL)
        if outflow_cut == true
            problem = include_cut_cell(problem, cut_cell_factor, 3)
            problem = include_cut_cell(problem, 0.001, 6)
            problem = include_cut_cell(problem, 0.45, 9)
            problem = include_cut_cell(problem, 0.1, 12)
            # additional cut cells
            problem = include_cut_cell(problem, 0.15, 15)
            problem = include_cut_cell(problem, 0.2, 18)
            problem = include_cut_cell(problem, 0.25, 21)
            problem = include_cut_cell(problem, 0.3, 24)
            problem = include_cut_cell(problem, 0.35, 27)
            problem = include_cut_cell(problem, 0.4, 30)
            problem = include_cut_cell(problem, 0.49, 33)
            problem = include_cut_cell(problem, 0.05, 36)
            problem = include_cut_cell(problem, 0.01, 39)
        else
            problem = include_cut_cell_only(problem, cut_cell_factor, 3)
        end
        RHS_mat, problem = DGsemidiscretization_DoD(problem, order, basis, "Upwind", do_stabilize = true, fix_eta = fix_eta, c = c)
        solution, u_exact = TMM(problem, RHS_mat)
        #errors[istep] = norm(solution[:, end] - u_exact[:, end])
         
        x_d = problem["x_d"]
        v = problem["v"]
        deg = problem["deg"]
        nodes = problem["nodes"]
        errors[istep] = 0
        for i in 1:problem["cellnumber"]
            basis1 = nodes(deg)
            for j in 1:deg+1
                # Inserting these nodes into the structure basis1, to adapt the field basis1.interpolation_matrix to the used nodes.
                basis1.nodes[j] = x_d[(deg+1)*(i-1) + j] 
            end
            errors[istep] += integrate((solution[(deg+1)*(i-1)+1:(deg+1)*i, end] - u_exact[(deg+1)*(i-1)+1:(deg+1)*i, end]).^2, basis1.weights.*(v[i+1]-v[i])/2)
        end
        errors[istep] = sqrt(errors[istep])
        
    end
    return stepsizes, errors
end

# Compare the error of the background method, the classic CFL-based choice for eta and the optimized choice, in every case with appropiate sharp CFL condition 
# (currently up to order 4, requires a new timestep method/the appropiate lambda for order 5) 
# JUST FOR CUT CELL FACTOR = 0.3! Other choices require lower CFL numbers for the background method, which makes running the convergence analysis in that case really slow
function convergence_error_comparison(basis, order; exprange = [3,8], sec_factor = 1.0, alignCFL = false, Tmax = 5)
    # LobattoLegendre denotes the first row of lambdas/CFL, GaussLegendre the second row of lambdas/CFL
    if basis == LobattoLegendre
        nodes_index = 1
    elseif basis == GaussLegendre
        nodes_index = 2
    end
    lambdas_opt = [1 0.8767 0.5398 0.3213 ; 1 0.7891 0.4416 0.2787]
    CFL_backgr = [1 1 0.44999 0.74672 ; 1 0.33481 0.2098 0.449994]
    CFL_DoD_classic = [1 .789668 0.164004 0.0976313 ; 1 0.323104 0.119105 0.0815261]
    CFL_opt = [1 0.78966 0.27625 0.358243 ; 1 0.32700 0.152292 0.26454]
    tmm = [expl_Euler, Heun, SSPRK3, SSPRK10_4]
    if alignCFL == true
        CFL_min = zeros(length(CFL_opt))
        for i in range(1,length(CFL_opt))
            CFL_min[i] = minimum([CFL_backgr[i], CFL_DoD_classic[i], CFL_opt[i]])
            CFL_backgr[i] = CFL_min[i]
            CFL_DoD_classic[i] = CFL_min[i]
            CFL_opt[i] = CFL_min[i]
        end
    end

    steps = zeros(exprange[2]-exprange[1]+1, 3)
    errors = zeros(exprange[2]-exprange[1]+1, 3)
    @time steps[:, 1], errors[:, 1] = test_convergence_background_method(basis, order-1; exprange = exprange, TMM = tmm[order], CFL = sec_factor*CFL_backgr[nodes_index, order], a=1, Tmax = Tmax)
    @time steps[:, 2], errors[:, 2] = test_convergence_dod_method(basis, order-1; exprange = exprange, TMM = tmm[order], CFL = sec_factor*CFL_DoD_classic[nodes_index, order], a=1, outflow_cut = true, fix_eta = false, cut_cell_factor = 0.3, Tmax = Tmax)
    @time steps[:, 3], errors[:, 3] = test_convergence_dod_method(basis, order-1; exprange = exprange, TMM = tmm[order], CFL = sec_factor*CFL_opt[nodes_index, order], a=1, outflow_cut = true, fix_eta = true, c = lambdas_opt[nodes_index, order], cut_cell_factor = 0.3, Tmax = Tmax)

    return steps, errors
end

#############################################################################
######################## Operator Norm calculation ##########################
#############################################################################

function calc_op_norm(setup, RHS_mat)
    x_d = setup["x_d"]
    M_global = setup["M"]
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
    return opnorm(sqrtM * RHS_mat * invsqrtM)
end

function op_norm_vs_cut_fac(deg, nodes, c_factors; cellnumber = 50, cutnumber = 3, with_cut_cell = false, do_stabilize = false, fix_eta = true, outflow_cut = true, CFL = 0.9, c = 1)
    opnorms = zeros(length(c_factors))
    for (ialpha,alpha) in enumerate(c_factors)
        problem = setup_problem_eq("sin", 0, 1, cellnumber, Tmax = 5, a = 1, CFL = CFL)
        if with_cut_cell
            if outflow_cut == true
                for cutn in cutnumber
                    problem = include_cut_cell(problem, alpha, cutn)
                    #problem = include_cut_cell(problem, alpha, 28)
                end
            else
                for cutn in cutnumber
                    problem = include_cut_cell_only(problem, alpha, cutn)
                end
            end
        end
        
        RHS_mat, problem = DGsemidiscretization_DoD(problem, deg, nodes, "Upwind", do_stabilize = do_stabilize, fix_eta = fix_eta, c = c)

        ###Option csv-file-input
        #RHS_mat = readdlm("./RHS_mat_legendre/M_$(ialpha).csv", ',', Float64)
        ###
        
        opnorms[ialpha] = calc_op_norm(problem, RHS_mat)
    end
    return opnorms
end


#############################################################################
############################# Further functions #############################
#############################################################################

# Full theoretical operator norm function from the eigenvalue analysis with Gerschgorin-circles for p = 0
function opnorm_theoretical(c, x)
    #c is the chosen parameter in 1-eta = abs(E)/(dt*a), if we choose dt = c*dx/a instead of the CFL-dependant timestep dt.
    # x is equal to the cut-cell factor alpha
    t1 = 1 + 1/c^2 + (1-x/c)^2/(1-x) + 1 + (1-x/c) + 1/c * abs((c-x)/((1-x)*c)-1/(c^2))
    t2 = 1/c^2*(1+x) + 1/((1-x)*c) * sqrt(x/(1-x)) + 1/c * abs((c-x)/((1-x)*c)-1/(c^2))
    t3 = 2/(1-x)^2 + 1/sqrt(1-x)+(1-x/c)/sqrt(1-x) + 1/c*sqrt(x/(1-x))
    return sqrt(maximum([t1, t2, t3]))
end

# Full theoretical operator norm function from the direct estimation for p = 0

function opnorm_theoretical_direct(c, x; optimized = true)
    #c is the chosen parameter in 1-eta = abs(E)/(dt*a), if we choose dt = c*dx/a instead of the CFL-dependant timestep dt.
    # x is equal to the cut-cell factor alpha
    
    fcut(x, xi1, xi2, ep1, ep2) = (1+1/xi2)/(c^2) - 4 + 2*x/(c^2)*(1 + 1/ep1 + 1/xi1)
    fm1(x, xi1, xi2, ep1, ep2) = x*(1+xi2)/(c^2)-2+2*(1-x/c)^2+2*(1-x)^2/ep2+2*(1-x)^2*xi1
    fp1(x, xi1, xi2, ep1, ep2) = 2*((1+ep1+ep2)/(1-x)-1)

    # target = [xi1, xi2, q1, q2]
    if optimized
        # old and possibly wrong, when fcut had an xi2 instead of xi1 in the end
        #target = [0.01, 3.0, 0.3422222222222222, 0.916060606060606]
        target = [0.6442424242424243, 2.4865656565656566, 0.463030303030303, 0.6744444444444444]
    else
        target = [1, 1, 1, 1]
    end
    
    xi1 = target[1]
    xi2 = target[2]
    ep1 = target[3] * x
    ep2 = target[4] * (1-x)

    t1 = fcut(x, xi1, xi2, ep1, ep2)
    t2 = fm1(x, xi1, xi2, ep1, ep2)
    t3 = fp1(x, xi1, xi2, ep1, ep2)


    return sqrt(4+maximum([t1, t2, t3]))
end

function opnorm_theoretical_p1_Lobatto(c, x; receive = "opnorm")
#c is the chosen parameter in 1-eta = abs(E)/(dt*a), if we choose dt = c*dx/a instead of the CFL-dependant timestep dt.
    # x is equal to the cut-cell factor alpha
    xi2= sqrt(5)-2    
    B = 2 + 2/xi2
    eta = 1 - x/c
    alpha = x

    # partial variable coefficients
    uim2_2 = 20-(6+2*xi2)
    uim1_1 = -B+6*(1+eta*alpha^2)^2 + 30*(eta*alpha)^2 + (1-eta*alpha*(2+alpha))^2 + 2*eta^2 + 4/(1-alpha)*(alpha*eta)^2
    uim1_2 = -B + 2 + 8*(eta*alpha*(2+alpha))^2 + (2*(1-eta)-eta*alpha)^2 + (eta*alpha)^2 + 16*((1+alpha)*eta)^2 + 4*((1+eta)*(1-alpha))^2 
            + (((1+alpha)/(1-alpha))*eta)^2 + 5*(1-eta) + eta*alpha
    ui_1 = -B + 2*eta + ((1-eta)/alpha)^2 + 4*((1-eta)/alpha)^3 + 9*eta^2*alpha
    ui_2 = -B + 2*eta + (6+4*eta)*((1-eta)/alpha)^2 + (12+4/(1-alpha))*((1-eta)^2/alpha) + 9*eta^2*alpha + ((1-eta)/alpha)^3 
            + 4*((1-eta)/(1-alpha))^2
    uip1_1 = -B + 2/((1-alpha)^2) + 3/((1-alpha)^3)
    uip1_2 = -B + (6+2*xi2)/(1-alpha) + 2/((1-alpha)^2) + 3/((1-alpha)^3)

    # prepare output
    sv_dict = Dict("uim2_2" => uim2_2, "uim1_1" => uim1_1, "uim1_2" => uim1_2, "ui_1" => ui_1, "ui_2" => ui_2,
                        "uip1_1" => uip1_1, "uip1_2" => uip1_2 )
    single_vals = [sv_dict["uim2_2"], sv_dict["uim1_1"], sv_dict["uim1_2"], sv_dict["ui_1"], sv_dict["ui_2"], sv_dict["uip1_1"], sv_dict["uip1_2"]]
    
    # output
    if receive == "opnorm"
        return sqrt(B+maximum(abs.(single_vals)))
    else
        return sv_dict[receive]
    end
end
function opnorm_theoretical_p1_Gauss(c, x; receive = "opnorm")
#c is the chosen parameter in 1-eta = abs(E)/(dt*a), if we choose dt = c*dx/a instead of the CFL-dependant timestep dt.
    # x is equal to the cut-cell factor alpha
    if c != 1.0
        println("p=1 Gauss Legendre just for c=1 implemented yet.")
    end
    Q = sqrt(3)
    B = 26 + 14*Q
    alpha = x
    eta = 1 - alpha


    # partial variable coefficients
    uimm_1 = (-8.5 + Q/2 + 8 + 4*Q 
        + 4+2*Q 
        + 7+3*Q + alpha^2*eta^2*(3*(2+(3+Q)*alpha))^2
        + 3+3*Q + eta^2*alpha^2*(3*alpha*(1+Q)+6*Q+9)^2
        + alpha*eta^2*(1+Q)^2
        + alpha*eta^2*(1+Q)^2)
    uimm_2 = (-(20+8*Q) + 8+4*Q
        + 4+2*Q
        + 5+3*Q + 3^2*alpha^2*eta^2*(2+(Q+3)*alpha)^2
        + 3+3*Q + alpha^2*eta^2*(4*Q-6 + 3*alpha*(1+Q))^2
        + alpha*eta^2*(1+Q)^2
        + alpha*eta^2*(1+Q)^2)    
    uim_1 = (-B + 6*alpha^2*(1-alpha)^2*(3*alpha^2+(6-2*Q)*alpha + 4-2*Q)+6*alpha*(1-alpha)^2 + 
        + 6*alpha^2*(1-alpha)*(7+Q)+(33-7*Q)*alpha*(1-alpha) +  16 - 2*Q
        + 7 + 3*Q + 1
        + 5+3*Q + 1
        + 6*alpha^2*(1-alpha)*(Q+4) + 2*alpha*(1-alpha)*(2*Q+9) + 4*(1-2*alpha)
        + 6*Q*alpha^3*eta^2 + 6*alpha^2*eta^2*(Q-1) +6*alpha*eta*(Q+1)+alpha*(6*alpha + Q)
        + 6*Q*alpha^3*eta^2 + 6*alpha^2*eta^2*(Q-1) +6*alpha*eta*(Q+1) + alpha*(Q*(1-2*alpha) + 8)
        + 6*alpha+(3-Q)
        + 2*alpha*(2*Q-3)+3*Q-5)
    uim_2 = (-B + 18*alpha^4*(1-alpha)^2 + 12*alpha^3*(1-alpha)^2(3+Q) + 12*alpha^2*(1-alpha)^2*(2+Q)
        + 6*alpha^2*(1-alpha)*(7-Q) + 6*alpha*(1-alpha)^2 + alpha*(1-alpha)*(30+6*Q) + 16+2*Q
        + 3+3*Q + 1
        + 3+Q + 1
        + 6*alpha^2*(1-alpha)*(Q+4) + 2*alpha*(1-alpha)*(2*Q+9) + 4*(1-2*alpha)
        + 4*Q*alpha^2 + 3*alpha + 3*alpha^2*(Q-1)
        + 2*alpha*(1-alpha)*(Q+3 + alpha*(4+Q) + 3*alpha^2*(1+Q))
        + 6*alpha + (Q+3) 
        + 2*alpha*(2*alpha-3)+Q-1)
    ui_1 =  (-B + 6*alpha*(1-alpha)^2 + 2*(4-Q) + 5/(1-alpha)
        + 1
        + 1
        + 6*Q*alpha^2*(1-alpha)^2 + 6*alpha*(1-alpha)^2*(Q-1) + 6*(1-alpha)*(Q+1) + (6*alpha+Q)
        + 4*Q*alpha + 3 + 3*alpha*(Q-1) 
        + 2*alpha*(abs(3*(1-alpha)^3-2)/(1-alpha))
        + 2/(1-alpha)*(2+Q)
        + (3+Q)/(1-alpha))
    ui_2 = (-B + 6*alpha*(1-alpha)^2 + 2*(4-Q) + 5/(1-alpha)
        + 1
        + 1
        + 6*Q*alpha^2*(1-alpha)^2 + 6*alpha*(1-alpha)^2*(Q-1) + 6*(1-alpha)*(Q+1) + Q*(1-2*alpha)+8
        + 2*(1-alpha)*(3+Q) + (4+Q) + 3*alpha*(1+Q)
        + 2*alpha*(abs(3*(1-alpha)^3-2)/(1-alpha))
        + (5+Q)/(1-alpha)
        + 2*Q/(1-alpha))
    uip_1 = (-B + (8.5-Q/2)/(1-alpha) + (2*Q+8)/(1-alpha)^2 
        + (6*alpha+(3-Q))/(1-alpha) 
        + (6*alpha+Q+3)/(1-alpha)
        + 2*alpha/(1-alpha)^2 *(2+Q) + alpha/(1-alpha)^2 *(5+Q) + 4*Q/(1-alpha)^2)
    uip_2 = (-B + (20+8*Q)/(1-alpha) + (8-2*Q)/(1-alpha)^2 
        + 2*alpha/(1-alpha)*(2*Q-3) + 1/(1-alpha)*(3*Q-5)
        + 2*alpha/(1-alpha)*(2*Q-3) + (Q-1)/(1-alpha)
        + alpha/(1-alpha)^2*(Q+3)
        + 2*alpha/(1-alpha)^2*Q
        + 4*Q/(1-alpha)^2)


    # prepare output
    sv_dict = Dict("uimm_1" => uimm_1, "uimm_2" => uimm_2, "uim_1" => uim_1, "uim_2" => uim_2, "ui_1" => ui_1, "ui_2" => ui_2,
                        "uip_1" => uip_1, "uip_2" => uip_2 )
    single_vals = [sv_dict["uimm_1"], sv_dict["uimm_2"], sv_dict["uim_1"], sv_dict["uim_2"], sv_dict["ui_1"], sv_dict["ui_2"], sv_dict["uip_1"], sv_dict["uip_2"]]
    
    # output
    if receive == "opnorm"
        return sqrt(B+maximum(abs.(single_vals)))
    else
        return sv_dict[receive]
    end
end
#############################################################################
############################ L^2-Norm over time #############################
#############################################################################

function l2_norm_over_time(basis, order , TMM; cellnumber = 50, CFL = 0.9, a=1, background = false, outflow_cut = true)
    problem = setup_problem_eq("sin", 0, 1, cellnumber, Tmax = 5, a = a, CFL = CFL)
    errors = zeros(length(problem["t_d"]))
    if background == false
        if outflow_cut == true
            problem = include_cut_cell(problem, 0.01, Integer(round(cellnumber/4)))
            problem = include_cut_cell(problem, 0.01, Integer(round(cellnumber/3)))
            problem = include_cut_cell(problem, 0.01, Integer(round(cellnumber/2)))
        else
            problem = include_cut_cell(problem, 0.01, Integer(round(cellnumber/4)))
            problem = include_cut_cell(problem, 0.01, Integer(round(cellnumber/3)))
            problem = include_cut_cell(problem, 0.01, Integer(round(cellnumber/2)))
        end
    end
    RHS_mat, problem = DGsemidiscretization_DoD(problem, order, basis, "Upwind", do_stabilize = true)
    solution, u_exact = TMM(problem, RHS_mat) 
    x_d = problem["x_d"]
    v = problem["v"]
    deg = problem["deg"]
    nodes = problem["nodes"]
    for (it, t) in enumerate(problem["t_d"])
        errors[it] = 0
        for i in 1:problem["cellnumber"]
            basis1 = nodes(deg)
            for j in 1:deg+1
                # Inserting these nodes into the structure basis1, to adapt the field basis1.interpolation_matrix to the used nodes.
                basis1.nodes[j] = x_d[(deg+1)*(i-1) + j]
            end
            errors[it] += integrate((solution[(deg+1)*(i-1)+1:(deg+1)*i, it]).^2, basis1.weights.*(v[i+1]-v[i])/2)
        end
        errors[it] = sqrt(errors[it])
    end
        
    return problem["t_d"], errors
end

#############################################################################
############# Determining long time stable CFL vs alpha #####################
#############################################################################

function calc_lts_CFL(basis, order , TMM; cellnumber = 50, a=1, background = false, fix_eta = false, outflow_cut = true, c = 1.0)
    #alphas = [0.01, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.425, 0.45, 0.475, 0.49]
    alphas = [0.01, 0.05, 0.075, 0.1, 0.1125, 0.125, 0.1375, 0.15,0.1625,0.175,0.1875,0.2,0.2125,0.225,0.2375,0.25,0.2625,0.275,0.2875,0.3,0.3125,0.325,0.3375,0.35,0.3625,0.375,0.3875,0.4,0.4125,0.425,0.4375,0.45, 0.475, 0.49]
    if background == true
        a_len_origin = length(alphas)
        alphas = [0.01]
    end
    CFL_vs_alpha = zeros(length(alphas))
    for (ialpha, alpha) in enumerate(alphas)
        println("Starting for alpha = $(alpha)")
        CFL_max = 2
        CFL_min = 0.001
        CFL_tol = 0.01

        CFL = CFL_max
        @time while CFL_max - CFL_min > CFL_tol*CFL
            CFL = (CFL_max+CFL_min)/2
            problem = setup_problem_eq("sin", 0, 1, cellnumber, Tmax = 50, a = 1, CFL = CFL)
            if background == false
                if outflow_cut == true
                    problem = include_cut_cell(problem, alpha, Integer(round(cellnumber/2))+1)
                else
                    problem = include_cut_cell_only(problem, alpha, Integer(round(cellnumber/2))+1)
                end
            end
            RHS_mat, problem = DGsemidiscretization_DoD(problem, order, basis, "Upwind", do_stabilize = true, fix_eta = fix_eta, c = c)
            solution, u_exact = TMM(problem, RHS_mat) 
            if (solution[:,end] .<= 1.005) == ones(size(solution[:, end]))
                CFL_min = CFL
            else
                CFL_max = CFL
            end
        end
        println("=> CFL = $(CFL)")
        CFL_vs_alpha[ialpha] = CFL
    end
    if background == true
        return range(0.01, 0.49, length = a_len_origin), ones(a_len_origin).*CFL_vs_alpha
    else
        return alphas, CFL_vs_alpha
    end
end

#############################################################################
############## Determining c that minimize operatornorm #####################
#############################################################################

function minimize_operatornorm(basis, deg; cellnumber = 50, a = 1, outflow_cut = true, verbose = true)
    alphas = range(0.001, 0.499, length = 51)
    cutnum = 26

    # CFL has no impact, but needs to be defined here
    CFL = 1.0
    
    # initializing
    best_norm = 0
    best_c = 0
    new_cs = zeros(25)
    ######################
    #first step
    ######################
    cs = range(0.01, 1, length = 51) 
    norm_matrix_vs_c = zeros(length(alphas), length(cs))
    # calculate all operatornorms
    if verbose == true
        println("Time for first 100 evaluations:")
        @time for (ic, c) in enumerate(cs)
            norm_matrix_vs_c[:, ic] = op_norm_vs_cut_fac(deg, basis, alphas; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, outflow_cut = outflow_cut, CFL = CFL, c = c)
        end
    else
        for (ic, c) in enumerate(cs)
            norm_matrix_vs_c[:, ic] = op_norm_vs_cut_fac(deg, basis, alphas; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, outflow_cut = outflow_cut, CFL = CFL, c = c)
        end
    end
    # calculate c for all alphas
    (best_norm, ind_best_c) = findmin(maximum(norm_matrix_vs_c, dims = 1))
    best_c = cs[ind_best_c[2]]
    if verbose == true
        println(L"First approximation result: $\|\|L\|\|$ = " * "$(best_norm), " * L"$c = $" * "$(best_c)")
    end
    ######################
    # iterate to more decimals
    ######################
    for exp in [-2, -3, -4]
        new_cs = range(best_c - 10.0^exp, best_c + 10.0^exp, length = 25)
        norm_matrix_vs_c = zeros(length(alphas), length(new_cs))
        if verbose == true
            println("Time for step $(-exp):")
            @time for (ic, c) in enumerate(new_cs)
                norm_matrix_vs_c[:, ic] = op_norm_vs_cut_fac(deg, basis, alphas; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, outflow_cut = outflow_cut, CFL = CFL, c = c)
            end
        else
            for (ic, c) in enumerate(new_cs)
            norm_matrix_vs_c[:, ic] = op_norm_vs_cut_fac(deg, basis, alphas; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, outflow_cut = outflow_cut, CFL = CFL, c = c)
            end
        end
       (best_norm, ind_best_c) = findmin(maximum(norm_matrix_vs_c, dims = 1))
        best_c = new_cs[ind_best_c[2]]
        if verbose == true
            println("Approximation result step $(-exp):" * L"$\|\|L\|\|$ = " * "$(best_norm), " * L"$c = $" * "$(best_c)")
        end
    end
    return best_norm, best_c
end


function minimize_operatornorm_cdep(basis, deg; cellnumber = 50, a = 1, outflow_cut = true, verbose = true)
    alphas = range(0.001, 0.499, length = 100)
    cutnum = 26

    # CFL has no impact, but needs to be defined here
    CFL = 1.0
    
    cs = range(0.01, 1, length = 100) 
    norm_matrix_vs_c = zeros(length(alphas), length(cs))
    # calculate all operatornorms
    if verbose == true
        println("Time for the 100 operatornorm evaluations:")
        @time for (ic, c) in enumerate(cs)
            norm_matrix_vs_c[:, ic] = op_norm_vs_cut_fac(deg, basis, alphas; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, outflow_cut = outflow_cut, CFL = CFL, c = c)
        end
    else
        for (ic, c) in enumerate(cs)
            norm_matrix_vs_c[:, ic] = op_norm_vs_cut_fac(deg, basis, alphas; with_cut_cell = true, do_stabilize = true, cutnumber = cutnum, fix_eta = true, outflow_cut = outflow_cut, CFL = CFL, c = c)
        end
    end
    
    (best_norm, ind_best_c) = findmin(norm_matrix_vs_c, dims = 2)
    
    cvals = zeros(length(cs))
    for k in 1:length(cs)
        cvals[k]=cs[ind_best_c[k][2]]
    end
    
    return alphas, best_norm, cvals
end






