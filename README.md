# TODO
- DOI
- references on paper and for this repo
- include instruction to obtain julia results
- include full instruction for DUNE results
- empty folders?

# The domain-of-dependence stabilization for cut-cell meshes is fully discretely stable

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](#TODO)

This repository contains information and code to reproduce the results
presented in the article
```bibtex
@online{
}
```

If you find these results useful, please cite the article mentioned above.
If you use the implementations provided here, please **also** cite this
repository as
```bibtex
@misc{
}
```

## Abstract

We present a fully discrete stability analysis of the
domain-of-dependence stabilization for hyperbolic problems. The method
aims to address issues caused by small cut cells by redistributing
mass around the neighborhood of a small cut cell at a semi-discrete
level. Our analysis is conducted for the linear advection model
problem in one spatial dimension. We demonstrate that fully discrete
stability can be achieved under a time step restriction that does not
depend on the arbitrarily small cells, using an operator norm
estimate. Additionally, this analysis offers a detailed understanding
of the stability mechanism and highlights some challenges associated
with higher-order polynomials. We also propose a way to mitigate these
issues to derive a feasible CFL-like condition.
The analytical findings, as well as the proposed
solution are verified numerically in
one- and two-dimensional simulations.


## Numerical experiments

In order to generate the results from this repository, you need to install [Julia](https://julialang.org).
We recommend using `juliaup`, as detailed in the official website [https://julialang.org](https://julialang.org).


The one dimensional results have been generated using Julia version 1.10.2, and we recommend installing the same.
Once you have installed Julia, you can clone this repository, enter this directory and start the executable
`julia` with the following steps

```shell
git clone https://github.com/louispetri/2025_dod_linear_stability.git
cd 2025_dod_linear_stability
julia --project=.
```

Then enter the following commands to generate the data for the one dimensional case and all the plots of the paper.

```julia
julia> import Pkg; Pkg.instantiate() # Does not need to be re-run the next time you enter the REPL
julia> include("--.jl") # What gets generated
# TODO
```

All the figures are now ready and available in the following locations:
1. 

The data for the two dimensional case were obtained using the DUNE framework in C++ and can be reproduced as follows:

Choose the desired polynomial degree. For this you can set the global variable order in cplusplus_code/dune-cutcell/src/simulationsscalartransport.cc.
Execute the file build.sh inside the directory cplusplus_code. It will download and build the necessary DUNE modules
and the TPMC package. Also make sure that GNU parallel (https://www.gnu.org/software/parallel/) is available in your path.
Move to the folder cplusplus_code/dune-cutcell/build. It should contain an executable src/simulationsscalartransport and
a configuration file scalar-transport.ini.
Inside the folder cplusplus_code/dune-cutcell/build you should execute the scripts cplusplus_code/dune-cutcell/utility/compute-operator-norms.sh
and cplusplus_code/dune-cutcell/compute-errors.sh.
Both scripts take two arguments, a file containing the desired offsets for the channel geometry (you can take the file located
at cplusplus_code/dune-cutcell/utility/offsets.txt) and the maximum number of jobs for GNU parallel. You can pass 0 to let
GNU parallel run as many jobs in parallel as possible.
The script cplusplus_code/dune-cutcell/utility/compute-operator-norms.sh will compute different operator norms. The outer loop
controls the different capacity factors that will be used. Make sure that the key output.operatorNorm in scalar-transport.ini
is set to true. The keys xresolution and channelTest.channelAngle, resp. channelTest.vfAngle control the number of background cells
and the angle of the channel. To reproduce operator norms from the paper, set xresolution and yresolution to 20. Also set channelTest.length to 1.0.
The final output will be a folder for each capacity factor containing a JSON file for each offset containing the computed operator norm.
The script cplusplus_code/dune-cutcell/utility/compute-errors.sh performs simulations for given CFL factors. The cfl factors can
be controlled by the outer loop inside the script. Keep in mind that a CFL factor of 1 / (2p + 1) is applied regardless of the
specified CFL factor. Make sure that the keys output.l2Norm and output.linfNorm in scalar-transport.ini are set to true. To
reproduce results from the paper, set xresolution and yresolution to 120 and channelTest.length to 6.0. The timestepping scheme
can be selected via the key timestepping.method. Available schemes are given by expliciteuler, heun, shu3, spprk104. The simulation time
can be set by timestepping.T. Set it to 5.0 to reproduce results from the paper.


## Authors

- Gunnar Birke (gunnar.birke@uni-muenster.de)
- [Christian Engwer](https://www.uni-muenster.de/FB10srvi/persdb/MM-member.php?id=724) (christian.engwer@uni-muenster.de)
- Louis Petri (lpetri01@uni-mainz.de)
- [Hendrik Ranocha](https://ranocha.de) (hendrik.ranocha@uni-mainz.de)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
