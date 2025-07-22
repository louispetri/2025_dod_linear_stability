#!/bin/bash

exec=$1
config=$2
start=$3
end=$4
step=$5
gpjobs=$6
mpinp=$7

# shellcheck disable=SC2086
# shellcheck disable=SC2207
resolutions=($(seq $start $step $end))
parallel -j "$gpjobs" mpirun -np "$mpinp" "$exec" -config "$config" -xresolution {} ::: ${resolutions[*]}
