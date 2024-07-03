#!/bin/bash

cp EoS/scripts/EoS_I_G2.py EoS/scripts/EoS.py

dm_om_arr=(0.25 0.35355339059 0.5 0.70710678118 1 1.41421356237 2 2.82842712475 4)
M_arr=(250 353.55339059 500 707.10678118 1000 1414.21356237 2000 2828.42712475 4000)



python3 src/plotting.py 1000 0
python3 src/plotting.py 1000 -1
python3 src/plotting.py 1000 1

for dm_om in "${dm_om_arr[@]}"
do
for M in "${M_arr[@]}"
do
	python3 src/plotting.py $M $dm_om
done
done
