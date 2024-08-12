#!/bin/bash

python3 src/plot_3D.py 1 max_mass
python3 src/plot_3D.py 2 max_mass
python3 src/plot_3D.py 3 max_mass
python3 src/plot_3D.py 1 k2_crit
python3 src/plot_3D.py 2 k2_crit
python3 src/plot_3D.py 3 k2_crit
python3 src/plot_3D.py 1 k2_max
python3 src/plot_3D.py 2 k2_max
python3 src/plot_3D.py 3 k2_max
python3 src/plot_3D.py 1 r_crit
python3 src/plot_3D.py 2 r_crit
python3 src/plot_3D.py 3 r_crit
python3 src/plot_3D.py 1 comp
python3 src/plot_3D.py 2 comp
python3 src/plot_3D.py 3 comp




python3 src/plotting.py 1000 0
python3 src/plotting.py 1000 -1
python3 src/plotting.py 1000 1

dm_om_arr=(0.25 0.35355339059 0.5 0.70710678118 1 1.41421356237 2 2.82842712475 4)
M_arr=(250 353.55339059 500 707.10678118 1000 1414.21356237 2000 2828.42712475 4000)

for dm_om in "${dm_om_arr[@]}"
do
for M in "${M_arr[@]}"
do
	python3 src/plotting.py $M $dm_om
done
done
