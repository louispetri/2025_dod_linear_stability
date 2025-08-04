
# Create necessary folders to store data

mkpath("./Plot_scripts/paper_plots_with_makie/CFL_lts_data")
mkpath("./Plot_scripts/paper_plots_with_makie/GLL_opnorm_quotient_data")
mkpath("./Plot_scripts/convergence_sharp_results/convergence_for_paper/Sec_fac_0.95_aligncfl")
mkpath("./Plot_scripts/convergence_sharp_results/convergence_for_paper/conv_work_precision_data")
mkpath("./Plot_scripts/Work_Precision_Diagramms/work_precision_data")


# 1.: Operator norm analysis for the extrapolation operator (Figure 2)

include("./Plot_scripts/opnorm_vs_extrapolation/compute_opnorm_extrapolation.jl") # Computation/Plotting

# 2.: Operator norm plots for DoD with nonoptimized parameters lambda_c (Figure 3/Figure 4)

include("./Plot_scripts/paper_plots_with_makie/paper_plots_makie.jl") # Computation/Plotting

# 3.: Optimized values for lambda_c (Table 1)

include("Plot_scripts/lambda_c_analysis/get_lambda_c_values.jl") # Computation

# 5. Operator norm quotient plots (Figure 5)

include("./Plot_scripts/paper_plots_with_makie/opnorm_quotient_makie.jl") # Computation

include("./Plot_scripts/paper_plots_with_makie/opnorm_quotient_makie_plot.jl") # Plotting

# 6. Partial operator norm plots, i.e., the involved block matrices (Figure 6)

include("./Plot_scripts/paper_plots_with_makie/Blockmatrix_details/blockmatrix_calc_and_eval.jl") # Computation and plotting
	
# 7. Sharp CFL restrictions for different alpha (Figure 7)

include("./Plot_scripts/paper_plots_with_makie/CFL_lts_makie.jl") # Computation

include("./Plot_scripts/paper_plots_with_makie/CFL_lts_makie_plots.jl") # Plotting

# 8. Convergence errors for DoD, applying different lambda_c (Figure 8)

include("./Plot_scripts/convergence_sharp_results/convergence_for_paper/convergence_get_data_aligncfl.jl") # Computation (requires time)

include("./Plot_scripts/convergence_sharp_results/convergence_for_paper/convergence_get_plots_makie_aligncfl.jl") # Plotting

# 9. Work precision diagrams for DoD, applying different lambda_c (Figure 9)

include("./Plot_scripts/convergence_sharp_results/convergence_for_paper/convergence_get_data_work_precision.jl") # Computation of errors (requires time)

include("./Plot_scripts/Work_Precision_Diagramms/work_precision_makie.jl") # Compute necessary data

include("./Plot_scripts/Work_Precision_Diagramms/work_precision_plot_makie.jl") # Plotting

# 10. Operator norms plots for the ramp setup with a 45° angle (Figure 11)
# Uses the data in the folder Julia_code/Plot_scripts/2d_plots_imported/Plotdaten.
# This data was obtained via DUNE and is just plotted in Julia.

include("./Plot_scripts/2d_plots_imported/opnorm_plots_2d_makie.jl")