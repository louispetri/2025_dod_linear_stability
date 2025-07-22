#!/bin/bash
set -x

offsets=$1
max_nodes=$2

# Every value will be tested as a cfl safety factor
# Note that a factor of 1 / (2p + 1) is always applied
for val in 1.0 1.1 1.2 1.3 1.4 1.5 ; do
	export cflSafetyFactor=$val
	sed_subs_cflSafetyFactor="s/cflSafetyFactor =.*/cflSafetyFactor = $cflSafetyFactor/g"
	sed -i -e "$sed_subs_cflSafetyFactor" scalar-transport.ini
	mkdir cfl=$cflSafetyFactor

	parallel -j "$max_nodes" -a "$offsets" 'sed_subs="s/offset =.*/offset = {}/g";'\
	'cp scalar-transport.ini scalar-transport-{}.ini;'\
	'sed -i -e "$sed_subs" scalar-transport-{}.ini;'\
	'mkdir offset={};'\
	'cd offset={};'\
	'../src/simulations-scalar-transport ../scalar-transport-{}.ini;'\
	'cd ..;'\
	'mv offset={}/ cfl=$cflSafetyFactor/;'\
	'rm scalar-transport-{}.ini'
done
