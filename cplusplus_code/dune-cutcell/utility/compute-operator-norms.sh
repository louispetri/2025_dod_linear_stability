#!/bin/bash
set -x

offsets=$1
max_nodes=$2

# Every value will be tested as a factor for the capacity
for val in 0.16 0.161 0.162 ; do
	export tau=$val
	sed_subs_tau="s/tau =.*/tau = $tau/g"
	sed -i -e "$sed_subs_tau" scalar-transport.ini
	mkdir tau=$tau

	parallel -j "$max_nodes" -a $offsets 'sed_subs="s/offset =.*/offset = {}/g";'\
	'cp scalar-transport.ini scalar-transport-{}.ini;'\
	'sed -i -e "$sed_subs" scalar-transport-{}.ini;'\
	'mkdir offset={};'\
	'cd offset={};'\
	'../src/simulations-scalar-transport -config ../scalar-transport-{}.ini;'\
	'mv operator-norm* operator-norm-{};'\
	'cd ..;'\
	'mv offset={}/operator-norm-{} tau=$tau;'\
	'rm -rf offset={};'\
	'rm scalar-transport-{}.ini'
done
