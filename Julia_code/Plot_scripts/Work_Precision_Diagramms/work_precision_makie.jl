include("work_precision.jl")
using  LaTeXStrings, DelimitedFiles


for basis in [LobattoLegendre, GaussLegendre]
    for order in [1,2,3,4]
        if basis == LobattoLegendre
            basissstr = "Lobatto"
        else
            basissstr = "Gauss"
        end
        steps_and_error =  work_precision_diagramm(basis, order, sec_factor = 0.95, alignCFL = false, relative_path_folder = "../convergence_sharp_results/convergence_for_paper/conv_work_precision_data")
        writedlm(joinpath(@__DIR__,"./work_precision_data/work_precision_data_$(basissstr)_order=$(order).txt"), steps_and_error)
    end
end
